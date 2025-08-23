# TransCisPredict: A Comprehensive Framework for Protein Level Prediction

This repository contains the complete analysis pipeline and data for our comprehensive framework to predict protein levels that incorporates both cis and trans genetic variants to facilitate conducting proteome-wide association studies (PWAS) on a biobank scale.

## Project Overview

Genetic effects on complex traits are mediated through protein pathways, making proteome-wide association studies (PWAS) increasingly important for understanding disease etiology. However, many studies lack proteomic data. We developed this comprehensive framework to predict protein levels and enable PWAS in datasets with genomic data but without measured protein levels.

### Key Features

- **Comprehensive Variant Selection**: Incorporates both cis and trans genetic variants across the genome
- **Multiple Statistical Methods**: Compares BayesR, SuSiE, LASSO, and elastic net approaches
- **Cross-Validation Optimization**: Selects optimal method for each protein individually
- **Large-Scale Application**: Designed for biobank-scale prediction in populations without proteomic data
- **Validated Framework**: Demonstrates strong correlation between predicted and measured protein associations

## Repository Structure

```
manuscript/
├── README.md                    # This file
├── data/                        # Data files, annotations, and weights
│   ├── bim_files_20250823_111518.tar.gz  # Variant annotation files (.bim format)
│   └── weights/                 # Final prediction weights (to be populated)
│       └── README.md           # Documentation for prediction weights
└── scripts/                     # Analysis pipeline scripts
    ├── README.md               # Detailed analysis pipeline documentation
    ├── step1_data_processing/  # Raw data processing and QC
    ├── step2_covariate_regression/  # Covariate adjustment
    ├── step3_ld_block_selection/    # Genomic region selection
    ├── step4_cross_validation/      # Method evaluation via CV
    ├── step5_cv_evaluation/         # CV performance assessment
    ├── step6_whole_sample_analysis/ # Final model training
    ├── step7_population_prediction/ # Proteome prediction
    ├── step8_pwas_analysis/         # Association testing
    └── utilities/                   # Common functions and methods
```

## Quick Start

1. **Setup Environment**: Install required R packages (see `scripts/README.md`)
2. **Configure Paths**: Modify placeholder paths in each script's configuration section
3. **Extract Data**: Unzip variant annotation files in `data/` directory  
4. **Run Pipeline**: Execute steps 1-8 sequentially following the analysis pipeline

## Key Components

### Scripts Directory
Contains the complete 8-step analysis pipeline with comprehensive documentation and usage examples. Each step is self-contained with clear input/output specifications.

### Data Directory  
Contains all data files including:
- Variant annotation files (.bim format) with genomic positions, rsIDs, and allele information
- Final trained protein prediction weights that can be applied to external genotype data (to be populated)

## Usage

See `scripts/README.md` for detailed pipeline documentation, usage instructions, and step-by-step execution guidelines.

## Applications

This framework enables:
- **Protein Expression Prediction**: Predict heritable component of protein levels using only genotype data
- **PWAS in Non-Proteomic Cohorts**: Conduct proteome-wide association studies without measured protein data
- **Disease Mechanism Discovery**: Uncover molecular mechanisms underlying complex diseases
- **Drug Target Identification**: Identify potential therapeutic targets through protein-trait associations

## Citation

Please cite the associated manuscript when using these analysis scripts and prediction weights for research purposes.

## License

This repository is intended for academic research and manuscript reproducibility. Please respect data use agreements for UK Biobank data.