# GS_Comparison_with_NCG
Comparing Genomic Prediction Methods Using Normalized Cumulative Gain

This repository provides R code to benchmark and compare various genomic prediction algorithms using both Pearson’s correlation coefficient (r) and Normalized Cumulative Gain (NCG) as evaluation metrics.

The script includes functions to:

* Apply multiple genomic prediction methods (linear and machine learning-based)

* Perform repeated cross-validation

* Evaluate and visualize model performance using both r and NCG

For more information on the use of Normalized Cumulative Gain in genomic prediction, please refer to our forthcoming publication:

Manuscript in preparation – link will be added upon publication.

## Data Format
The required genotype and phenotype data should be provided in the .raw format, as used by [PLINK](https://www.cog-genomics.org/plink/1.9/formats).

## Example Datasets
The datasets analyzed in our study are available as compressed files in the [Datasets](https://github.com/FelixHeinrich/GS_Comparison_with_NCG/tree/main/Datasets) folder.

## License

This project is licensed under the **GPL-3.0 License** - see [LICENSE](LICENSE) for more information.

## Citation

Heinrich, F.; Lange, T.M.; Ramzan, F.; Gültas, M.; Schmitt, A.O. Normalized cumulative gain as an alternative evaluation measure for genomic selection models. Genet Sel Evol 2025, 57, 70. [https://doi.org/10.1186/s12711-025-01022-9](https://doi.org/10.1186/s12711-025-01022-9)
