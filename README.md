# TransCisPredict: A Comprehensive Framework for Protein Level Prediction

This repository contains the complete analysis pipeline for the comprehensive framework to predict protein levels that incorporates both _cis_- and _trans_- variants to facilitate conducting proteome-wide association studies (PWAS) on a biobank scale. The weights that have been generated from UK Biobank are also available through [Synapse ID: 69052240](https://www.synapse.org/Synapse:syn69052240/wiki/634058).

## Repository Structure

```
TransCisPredict/
├── README.md                    # This file
└── scripts/                     # Analysis pipeline scripts
    ├── README.md               # Detailed analysis pipeline documentation
    ├── step1_data_processing/  # Raw data processing and QC
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

a. **Setup Environment**: Install required R packages (see `scripts/README.md`)
b. **Configure Paths**: Modify placeholder paths in each script's configuration section
c. **Obtain Prediction Weights**: Either download pre-computed weights from [Synapse ID: 69052240](https://www.synapse.org/Synapse:syn69052240/wiki/634058), or generate custom weights by running steps 1-6 of the analysis pipeline
d. **Predict Protein Levels**: Use the obtained weights to predict protein expression levels in your target sample (see `step7_population_prediction/`)
e. **Perform PWAS**: Conduct proteome-wide association studies using the predicted protein levels (see `step8_pwas_analysis/`)

## Key Components

### Scripts Directory
Contains the complete analysis pipeline with comprehensive documentation and usage examples. Each step is self-contained with clear input/output specifications.

## Applications

This framework enables:
- **Weight Estimation**: Generate weight estimation between protein expression level and genetic vriants from the reference sample
- **Protein Expression Prediction**: Predict heritable component of protein levels using only genotype data from the target sample
- **PWAS in Non-Proteomic Cohorts**: Conduct PWAS

## Citation

Please cite the associated manuscript when using these analysis scripts and prediction weights for research purposes:
> Dong R, Lamb D, Wang GT, DeWan AT, Leal SM. Leveraging cis- and trans-variants to improve plasma protein level prediction for proteome-wide association studies (2025). Submitted.

## License

This repository is intended for academic research and manuscript reproducibility. Please respect data use agreements for UK Biobank data.
