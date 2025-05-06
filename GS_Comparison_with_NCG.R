#Loading required libraries####
library(rrBLUP)
library(data.table)
library(ggplot2)
library(doParallel)
library(GROAN)
library(BGLR)
library(ranger)
library(gbm)
library(xgboost)
library(dplyr)
library(tidyr)
library(ggpubr)
library(dunn.test)
library(multcompView)
library(ggtext)

#Perform genomic prediction with cross-validation####

#' Analyse a dataset with different algorithms for genomic prediction and write the predictions to a file
#'
#' @param rawPath Path to a PLINK .raw file that contains the genotype and phenotype information
#' @param outputPath Path to file where the results should be saved
#' @param repeatNumber Number of times the cross-validation should be repeated
#' @param foldNumber Number of folds for cross-validation
#'
performGenomicPredictionsWithDifferentAlgorithms = function(rawPath, outputPath, repeatNumber, foldNumber){
  algoVec = c("BLUP","BRR", "BA","BB","BL","BC", "RF", "GBM", "XGBoost") # Names of the algorithms (for further algorithms add their names here)
  #set.seed(42) #Uncomment to reproduce the published results
  rawData = fread(rawPath, data.table = F, header = T)
  snpData = subset(rawData, select = -c(FID,IID,PAT,MAT,SEX))
  rownames(snpData) = rawData$IID
  snpCount = ncol(snpData)-1
  sampleCount = nrow(snpData)
  sampleIndices = 1:sampleCount
  foldsIndices = cut(sampleIndices, breaks = foldNumber, labels = F)
  
  resultsDF = as.data.frame(matrix(nrow = sampleCount*repeatNumber*length(algoVec), ncol = 6))
  colnames(resultsDF) = c("ID","Algo","Repeat","Fold","Pheno","Pred")
  
  #Running genomic prediction with repeated cross-validation
  entryCounter = 1
  for(rep in 1:repeatNumber){
    cat("Repetition: ", rep, "\n")
    repeatIndices = sample(sampleIndices) #For each repetition shuffle the individuals
    snpData = snpData[repeatIndices,]
    for(fold in 1:foldNumber){
      cat("Fold: ", fold, "\n")
      testIndices = repeatIndices[which(foldsIndices == fold, arr.ind = T)]
      trainIndices = setdiff(repeatIndices, testIndices)
      trainData = snpData[trainIndices,]
      testData = snpData[testIndices,]
      for(algo in algoVec){
        cat("Algo: ", algo, "\n")
        if(algo == "BLUP"){
          train_BLUP = mixed.solve(y = trainData[,1], Z = subset(trainData, select =  -get("PHENOTYPE")), K = NULL, SE = FALSE, return.Hinv = FALSE)
          train_BLUP.effects = as.matrix(train_BLUP$u)
          test_genoTimesEffects = as.matrix(subset(testData, select =  -get("PHENOTYPE"))) %*% train_BLUP.effects
          predValues = test_genoTimesEffects + c(train_BLUP$beta)
        }
        if(algo == "BRR"){
          fm = BGLR(y = trainData[,1], ETA = list(list(X = subset(trainData, select =  -get("PHENOTYPE")), model = "BRR")), nIter = 500, burnIn = 100)
          predValues = fm$mu + as.vector(as.matrix(subset(testData, select =  -get("PHENOTYPE"))) %*% fm$ETA[[1]]$b)
        }
        if(algo == "BA"){
          fm = BGLR(y = trainData[,1], ETA = list(list(X = subset(trainData, select =  -get("PHENOTYPE")), model = "BayesA")), nIter = 500, burnIn = 100)
          predValues = fm$mu + as.vector(as.matrix(subset(testData, select =  -get("PHENOTYPE"))) %*% fm$ETA[[1]]$b)
        }
        if(algo == "BB"){
          fm = BGLR(y = trainData[,1], ETA = list(list(X = subset(trainData, select =  -get("PHENOTYPE")), model = "BayesB")), nIter = 500, burnIn = 100)
          predValues = fm$mu + as.vector(as.matrix(subset(testData, select =  -get("PHENOTYPE"))) %*% fm$ETA[[1]]$b)
        }
        if(algo == "BL"){
          fm = BGLR(y = trainData[,1], ETA = list(list(X = subset(trainData, select =  -get("PHENOTYPE")), model = "BL")), nIter = 500, burnIn = 100)
          predValues = fm$mu + as.vector(as.matrix(subset(testData, select =  -get("PHENOTYPE"))) %*% fm$ETA[[1]]$b)
        }
        if(algo == "BC"){
          fm = BGLR(y = trainData[,1], ETA = list(list(X = subset(trainData, select =  -get("PHENOTYPE")), model = "BayesC")), nIter = 500, burnIn = 100)
          predValues = fm$mu + as.vector(as.matrix(subset(testData, select =  -get("PHENOTYPE"))) %*% fm$ETA[[1]]$b)
        }
        if(algo == "RF"){
          forest = ranger(y = trainData[,1], x = subset(trainData, select = -get("PHENOTYPE")), num.threads = 20)
          predValues = predict(forest, testData)$predictions
        }
        if(algo == "GBM"){
          gbm_model = gbm.fit(y = trainData[,1], x = subset(trainData, select = -get("PHENOTYPE")), distribution = "gaussian")
          predValues = predict(gbm_model, testData)
        }
        if(algo == "XGBoost"){
          xgb_train = xgb.DMatrix(data = as.matrix(trainData[,-1]), label = trainData[,1])
          xgboost_model = xgboost(data = xgb_train, nround = 500, early_stopping_rounds=100, print_every_n=0, nthread=20, verbose = 0)
          predValues = predict(xgboost_model, xgb.DMatrix(data = as.matrix(testData[,-1])))
        }
        #
        # Further algorithms can be added here with additional if functions
        #
        trueValues = testData$PHENOTYPE
        for(i in 1:nrow(testData)){
          resultsDF[entryCounter,] = c(rownames(testData)[i], algo, rep, fold, trueValues[i], predValues[i])
          entryCounter = entryCounter+1
        }
      }
    }
  }
  write.table(resultsDF, outputPath)
}

#Apply function on the different datasets
performGenomicPredictionsWithDifferentAlgorithms("Datasets/Goats.raw", "Results/Goats.results", repeatNumber = 10, foldNumber = 5)
performGenomicPredictionsWithDifferentAlgorithms("Datasets/Chicken.raw", "Results/Chicken.results", repeatNumber = 10, foldNumber = 5)
performGenomicPredictionsWithDifferentAlgorithms("Datasets/Wheat.raw", "Results/Wheat.results", repeatNumber = 10, foldNumber = 5)
performGenomicPredictionsWithDifferentAlgorithms("Datasets/Rice.raw", "Results/Rice.results", repeatNumber = 10, foldNumber = 5)

#Visualize the results####

#' Visualize the Pearson's correlation coefficient values of the algorithms using boxplot with non-parametric tests for significant differences
#'
#' @param dataset Dataframe containing the correlation coefficients for different algorithms across all folds
#' @param model_colors Colors that should be used for the different algorithms
#' @param thresh p-value threshold to mark significant differences between the algorithms
#' @param titleText Title for the plot
#'
createMeanComparisonBoxplot_NonParametric = function(dataset, model_colors, thresh, titleText){
  D = dataset
  D$Algo = as.factor(D$Algo)
  # Perform Kruskal Wallis and Dunn
  kruskal.res1 = kruskal.test(Pearson ~ Algo, data = D)
  dunn.res = dunn.test(x = D$Pearson, g = D$Algo, list = T, method = "Holm")
  
  # Extract adjusted p-values and format into a data frame
  pval_df <- data.frame(Comparison = dunn.res$comparisons, P.value = dunn.res$P.adjusted)
  # Split "Algo1 - Algo2" into separate columns
  pval_df <- separate(pval_df, Comparison, into = c("Algo1", "Algo2"), sep = " - ")
  all_algos <- unique(c(pval_df$Algo1, pval_df$Algo2))
  # Convert to a distance matrix format for multcompLetters()
  p_matrix <- matrix(NA, nrow = length(all_algos), ncol = length(all_algos),
                     dimnames = list(all_algos, all_algos))
  for (i in seq_len(nrow(pval_df))) {
    algo1 <- pval_df$Algo1[i]
    algo2 <- pval_df$Algo2[i]
    p_matrix[algo1, algo2] <- pval_df$P.value[i]
    p_matrix[algo2, algo1] <- pval_df$P.value[i]
  }
  diag(p_matrix) <- 1
  # Generate compact letter display (CLD)
  cld <- multcompLetters(p_matrix, threshold = thresh)
  
  # Extract only the grouping letters and convert to a data frame
  groups_df <- data.frame(
    Algo = names(cld$Letters),  # Algorithm names
    Letter = cld$Letters        # Corresponding letter groups
  )
  
  # Merge CLD letters with dataset
  D <- D %>%
    left_join(groups_df, by = "Algo")
  
  # Plot with letters above boxes
  ggplot(D, aes(x = Algo, y = Pearson, fill = Algo)) +
    geom_boxplot() +
    scale_fill_manual(values = model_colors) +
    geom_text(aes(label = Letter, y = max(Pearson) * 1.05), size = 5) +  # Add letters above boxes
    labs(
      title = titleText,
      y = "Pearson's correlation coefficient",
      fill = "Algorithm",
      x = element_blank()
    ) +
    theme_minimal() +
    theme(text = element_text(size = 24), axis.text.x = element_blank(), axis.ticks.x = element_blank())
}

#' Visualize the mean NCG values of the algorithms across all selection sizes
#'
#' @param dataset Dataframe containing the NCG values for different algorithms across all folds and selection sizes
#' @param model_colors Colors that should be used for the different algorithms
#' @param titleText Title for the plot
#'
createMeanNCGPlot = function(dataset, model_colors, titleText){
  # 1. Calculate mean performance for each model
  mean_trends <- dataset %>%
    group_by(Algo, k) %>%
    summarize(Mean_nCG = mean(nCG), .groups = "drop")
  
  
  # 2. Find the best-performing model based on smoothed trends
  best_mean_model <- mean_trends %>%
    group_by(k) %>%
    summarize(
      Best_Algorithm = Algo[which.max(Mean_nCG)],
      Best_nCG = max(Mean_nCG)
    ) %>%
    ungroup()
  
  # 3. Create segments where the best-performing model remains constant
  best_segments <- best_mean_model %>%
    mutate(Segment = cumsum(Best_Algorithm != lag(Best_Algorithm, default = first(Best_Algorithm)))) %>%
    group_by(Segment) %>%
    summarize(
      Start = min(k),
      End = max(k),
      Best_Algorithm = first(Best_Algorithm)
    )
  
  # Adjust single-value segments: add small offsets
  epsilon <- 0.5  # Small width adjustment
  best_segments <- best_segments %>%
    mutate(
      xmin = ifelse(Start == End, Start - epsilon, Start),
      xmax = ifelse(Start == End, End + epsilon, End)
    )
  minMean = min(mean_trends$Mean_nCG)
  maxMean = max(mean_trends$Mean_nCG)
  
  # 4. Plot background segments with smoothed trend lines
  meanPlot = ggplot() +
    # Background rectangles for best model segments
    geom_rect(
      data = best_segments,
      aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = (maxMean-minMean)*0.1+minMean, fill = Best_Algorithm),
      alpha = 0.3
    ) +
    # Smoothed trend lines
    geom_line(
      data = mean_trends,
      aes(x = k, y = Mean_nCG, color = Algo),
      size = 1
    ) +
    # Customize colors and theme
    scale_fill_manual(values = model_colors) +
    scale_color_manual(values = model_colors) +
    labs(
      title = titleText,
      x = "Number of selected individuals (k)",
      y = "*NCG* - Mean",
      fill = "Best Algorithm",
      color = "Algorithm"
    ) +
    theme_minimal() + theme(text = element_text(size = 24), plot.title = element_markdown(), axis.title.y = element_markdown())
}


#' Calculate performance values and visualize them
#'
#' @param resultsPath Path to file where the results have been saved
#' @param species Name of species to be used for plot titles and filenames
#' @param trait Name of trait to be used for plot titles and filenames
#'
compareAndVisualizeDifferentAlgorithms = function(resultsPath, species, trait){
  resultDF = read.table(resultsPath, header = T)
  
  listNCG = list()
  listMeasures = list()
  entryCounter = 1
  for(algo in unique(resultDF$Algo)){
    for(rep in unique(resultDF$Repeat)){
      for(fold in unique(resultDF$Fold)){ # Calculate NCG and other measures for each fold of the repeated cross-validation
        currentData = subset(resultDF, resultDF$Algo == algo & resultDF$Repeat == rep & resultDF$Fold == fold)
        if(nrow(currentData) > 0){
          phenoOrderedByTruth = currentData$Pheno[order(currentData$Pheno, decreasing = T)]
          phenoOrderedByPred = currentData$Pheno[order(currentData$Pred, decreasing = T)]
          
          cumGain = cumsum(phenoOrderedByPred)
          maxCumGain = cumsum(phenoOrderedByTruth)
          nCG = cumGain / maxCumGain
          listNCG[[entryCounter]] = data.frame(Algo = rep(algo, nrow(currentData)), Rep = rep(rep, nrow(currentData)), Fold = rep(fold, nrow(currentData)), k = 1:nrow(currentData), nCG = nCG)
          
          measures = measurePredictionPerformance(currentData$Pheno, currentData$Pred)
          listMeasures[[entryCounter]] = data.frame(Algo = algo, Rep = rep, Fold = fold, Pearson = measures["pearson"], RMSE = measures["rmse"], CoeffDet = measures["coeff_det"],
                                                    NDCG10 = measures["ndcg10"], NDCG20 = measures["ndcg20"])
          entryCounter = entryCounter + 1
        }
      }
    }
  }
  ncgDF = do.call(rbind, listNCG) #NCG values for the different algorithms over all folds and selection sizes
  measuresDF = do.call(rbind, listMeasures) #Pearson's correlation coefficient and other measures for the different algorithms over all folds
  
  model_colors <- c(
    "BRR" = "red",
    "BA" = "blue",
    "BB" = "green",
    "BL" = "lightblue",
    "BC" = "darkgreen",
    "BLUP" = "orange",
    "RF" = "black",
    "GBM" = "yellow",
    "XGBoost" = "brown"
    # If you have added further algorithms you need to assign them colors here
  )
  
  #Create plots to compare the algorithms based on Pearson correlation coefficient or NCG
  pearsonComparisonPlot = createMeanComparisonBoxplot_NonParametric(measuresDF, model_colors, 0.05, 
                                            paste0("Performance of algorithms based on Pearson's correlation coefficient \u2013 ", species, " (", trait, ")"))
  ncgComparisonPlot = createMeanNCGPlot(ncgDF, model_colors, 
                                              paste0("Performance of algorithms based on *NCG* with band highlighting best algorithm \u2013 ", species, " (", trait, ")"))
  ggsave(paste0(species, "_", trait, "_Pearson_Comparison.png"), pearsonComparisonPlot, device = "png", width=1980/100, height=1280/100, dpi = 100)
  ggsave(paste0(species, "_", trait, "_NCG_Comparison.png"), ncgComparisonPlot, device = "png", width=1980/100, height=1280/100, dpi = 100)
}

compareAndVisualizeDifferentAlgorithms("Results/Murciano_Granadina_Goats_MY210_NoMissing.results", "goat", "milk yield")
compareAndVisualizeDifferentAlgorithms("Results/EW_36.results", "chicken", "eggweight")
compareAndVisualizeDifferentAlgorithms("Results/Wheat_Scott_Yield_17_18.results", "wheat", "yield")
compareAndVisualizeDifferentAlgorithms("Results/Rice_GY.results", "rice", "grain yield")
