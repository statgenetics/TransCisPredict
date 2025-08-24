# TransCisPredict: A Comprehensive Framework for Protein Level Prediction

This repository contains the complete analysis pipeline and trained weights for our comprehensive framework to predict protein levels that incorporates both cis and trans genetic variants to facilitate conducting proteome-wide association studies (PWAS) on a biobank scale. Due to size constraints, the complete datasets are available through Synapse (see Data Availability section below).

## Repository Structure

```
TransCisPredict/
├── README.md                    # This file
├── data/                        # Data availability information
│   └── README.md               # Instructions for accessing datasets via Synapse
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
3. **Access Data**: Follow instructions in `data/README.md` to download datasets from Synapse
4. **Run Pipeline**: Execute steps 1-8 sequentially following the analysis pipeline

## Key Components

### Scripts Directory
Contains the complete 8-step analysis pipeline with comprehensive documentation and usage examples. Each step is self-contained with clear input/output specifications.

### Data Directory  
Contains information for accessing datasets including:
- Instructions for downloading complete datasets from Synapse
- Dataset descriptions and file formats
- Access requirements and data use agreements

## Usage

See `scripts/README.md` for detailed pipeline documentation, usage instructions, and step-by-step execution guidelines.

## Applications

This framework enables:
- **Protein Expression Prediction**: Predict heritable component of protein levels using only genotype data
- **PWAS in Non-Proteomic Cohorts**: Conduct proteome-wide association studies without measured protein data
- **Disease Mechanism Discovery**: Uncover molecular mechanisms underlying complex diseases
- **Drug Target Identification**: Identify potential therapeutic targets through protein-trait associations

## Data Availability

Due to the large size of the datasets, all data files are hosted on Synapse. See `data/README.md` for detailed instructions on accessing:
- Variant annotation files (.bim format) with genomic positions, rsIDs, and allele information
- Final trained protein prediction weights for external genotype data application

## Citation

Please cite the associated manuscript when using these analysis scripts and prediction weights for research purposes.

## License

This repository is intended for academic research and manuscript reproducibility. Please respect data use agreements for UK Biobank data.