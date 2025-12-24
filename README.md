# TransCisPredict: A Comprehensive Framework for Protein Level Prediction

This repository contains the complete analysis pipeline for the comprehensive framework to predict protein levels that incorporates both _cis_- and _trans_- variants to facilitate conducting proteome-wide association studies (PWAS) on a biobank scale. The weights that have been generated from UK Biobank are also available through [Synapse ID: 69052240](https://www.synapse.org/Synapse:syn69052240/wiki/634058).

## Repository Structure

```
TransCisPredict/
├── README.md                    # This file
├── readme_scripts.md            # Detailed analysis pipeline documentation
├── step1_data_processing/       # Raw data processing and QC
├── step2_covariate_regression/  # Covariate adjustment
├── step3_ld_block_selection/    # Genomic region selection
├── step4_cross_validation/      # CV using BayesR, SuSiE, LASSO, and Elastic Net
├── step5_cv_evaluation/         # CV performance evaluation
├── step6_whole_sample_analysis/ # Weight Estimation using the "optimal" method
├── step7_population_prediction/ # Protein expression level prediction in the target sample
├── step8_pwas_analysis/         # Proteome-wide association analyses (PWAS)
└── utilities/                   # Commonly used functions
```

## Quick Start

a. **Setup Environment**: Install required R packages (see `readme_scripts.md`)

b. **Configure Paths**: Modify placeholder paths in each script's configuration section

c. **Obtain Prediction Weights**: Either download pre-computed weights from [Synapse ID: 69052240](https://www.synapse.org/Synapse:syn69052240/wiki/634058), or generate custom weights by running steps 1-6 of the analysis pipeline

d. **Predict Protein Levels**: Use the obtained weights to predict protein expression levels in your target sample (see step7 scripts)

e. **Perform PWAS**: Conduct proteome-wide association studies using the predicted protein levels (see step8 script)

## Analysis Pipeline
The complete analysis pipeline with comprehensive documentation and usage examples is available in `readme_scripts.md`. Each step is self-contained with clear input/output specifications.

## Applications

This framework enables:
- **Weight Estimation**: Generate weight estimation between protein expression level and genetic vriants from the reference sample
- **Protein Expression Prediction**: Predict heritable component of protein levels using only genotype data from the target sample
- **PWAS in Target Cohort**: Conduct PWAS to identify associations with complex traits

## Citation

Please cite the associated manuscript when using these analysis scripts and prediction weights for research purposes:
> Dong R, Lamb D, Wang GT, DeWan AT, Leal SM. Leveraging cis- and trans-variants to improve plasma protein level prediction for proteome-wide association studies (2025). Submitted.

## License

This repository is intended for academic research and manuscript reproducibility. Please respect data usage agreements for UK Biobank data [FIXME].
