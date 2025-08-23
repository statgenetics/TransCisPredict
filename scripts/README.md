# TransCisPredict Analysis Pipeline

This directory contains the complete analysis pipeline for the TransCisPredict framework - a comprehensive approach to protein expression imputation that incorporates both cis and trans genetic variants for proteome-wide association studies (PWAS) using UK Biobank Pharma Proteomics Project data. The pipeline is organized into 8 sequential steps, each contained in its own directory with self-contained R scripts.

## Pipeline Overview

### Step 1: Data Processing
**Directory**: `step1_data_processing/`
**Script**: `process_olink_data.R`

Processes raw OLINK proteomics data, performs quality control, and creates individual protein expression files.

### Step 2: Covariate Regression  
**Directory**: `step2_covariate_regression/`
**Script**: `regress_npx_covariates_residuals.R`

Removes covariate effects (age, sex, BMI, genetic PCs) from protein expression to generate residuals for genetic analysis.

### Step 3: LD Block Selection
**Directory**: `step3_ld_block_selection/`
**Script**: `ld_block_selection_by_fdr.R`

Identifies genomic regions with significant genetic signal using FDR-based selection to focus downstream analysis.

### Step 4: Cross Validation
**Directory**: `step4_cross_validation/`  
**Script**: `pqtl_analysis_with_cross_validation.R`

Performs cross-validation to evaluate multiple statistical methods (BayesR, LASSO, Elastic Net, SuSiE) for protein prediction using both cis and trans genetic variants.

### Step 5: CV Evaluation
**Directory**: `step5_cv_evaluation/`
**Scripts**: 
- `evaluate_cv_performance.R` - Evaluates cross-validation performance
- `identify_best_method.R` - Selects optimal method per protein  
- `summarize_all_methods.R` - Summarizes performance across all methods

Evaluates cross-validation results and identifies the best-performing statistical method for each protein.

### Step 6: Whole Sample Analysis
**Directory**: `step6_whole_sample_analysis/`
**Script**: `pqtl_analysis_whole_sample.R`

Trains final prediction models using the complete dataset and optimal method identified for each protein. Generates the final cis+trans prediction weights for proteins that meet performance thresholds.

### Step 7: Population Prediction
**Directory**: `step7_population_prediction/`
**Scripts**:
- `predict_npx_all_ukbb_white_europeans.R` - Predicts protein levels for target population
- `combine_npx_files.R` - Combines individual protein predictions into final matrix

Applies trained cis+trans models to predict protein levels for UK Biobank individuals without measured proteomic data.

### Step 8: PWAS Analysis  
**Directory**: `step8_pwas_analysis/`
**Script**: `250304_gee_pwas_analysis.R`

Performs proteome-wide association studies using generalized estimating equations to test associations between predicted proteins and complex traits.

## Utilities

**Directory**: `utilities/`
**Scripts**:
- `pqtl_functions.R` - Core utility functions for data processing and genetic prediction
- `pqtl_weights.R` - Statistical method implementations (BayesR, LASSO, Elastic Net, SuSiE) for both cis and trans variant incorporation
- `timing_function.R` - Runtime calculation utilities

Contains common functions and method implementations used throughout the TransCisPredict pipeline.

## Data Files

**Directory**: `../data/`
**Contents**:
- `weights/` - Final prediction weights for each protein from Step 6
  - `{protein_name}_final_weights.csv` - Weight files for individual proteins
  - `variants.bim` - Variant annotation file mapping variant IDs to genomic positions and rsIDs
- `bim_files_20250823_111518.tar.gz` - Compressed variant annotation files

The data directory contains trained protein prediction weights and variant annotations that can be applied to external genotype data for protein expression prediction.

---

## Usage Instructions

### Running the Pipeline

1. **Sequential Execution**: Run steps 1-8 in order, as each step depends on outputs from previous steps.

2. **Path Configuration**: Modify placeholder paths in each script's configuration section for your environment:
   - Input data directories
   - Output directories  
   - Genotype data locations
   - Reference files

3. **Command Line Usage**: Most scripts accept command line arguments. See individual script headers for specific usage examples.

4. **Example Workflow**:
   ```bash
   # Step 1: Process raw data
   Rscript step1_data_processing/process_olink_data.R
   
   # Step 2: Generate residuals (example for APOE protein)
   Rscript step2_covariate_regression/regress_npx_covariates_residuals.R APOE
   
   # Step 4: Cross-validation (example for APOE)
   Rscript step4_cross_validation/pqtl_analysis_with_cross_validation.R APOE input_dir output_dir 0.2
   
   # Continue through remaining steps...
   ```

### Requirements

- **R packages**: tidyverse, data.table, plink2R, pecotmr, geepack, broom.mixed, rsample
- **Computing resources**: Memory requirements vary by step

This pipeline provides a complete framework for protein expression imputation and proteome-wide association studies that incorporates both cis and trans genetic variants for improved prediction performance.