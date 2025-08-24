# TransCisPredict Analysis Pipeline

This directory contains the complete analysis pipeline for the TransCisPredict framework - a comprehensive approach to protein expression imputation that incorporates both cis and trans genetic variants for proteome-wide association studies (PWAS) using UK Biobank Pharma Proteomics Project data. The pipeline consists of 8 sequential steps, each implemented as a standalone R script with standardized configuration sections.

## Pipeline Overview and Usage

### Requirements and Setup

**Required R packages**: 
- [tidyverse](https://cran.r-project.org/package=tidyverse) - Data manipulation and visualization
- [data.table](https://cran.r-project.org/package=data.table) - High-performance data processing  
- [plink2R](https://cran.r-project.org/package=plink2R) - Genotype data handling
- [pecotmr](https://github.com/cumc/pecotmr) - Statistical methods implementation
- [rsample](https://cran.r-project.org/package=rsample) - Cross-validation and resampling
- [broom](https://cran.r-project.org/package=broom) - Model output formatting
- [broom.mixed](https://cran.r-project.org/package=broom.mixed) - Mixed model output formatting
- [janitor](https://cran.r-project.org/package=janitor) - Data cleaning utilities
- [geepack](https://cran.r-project.org/package=geepack) - Generalized estimating equations
- [igraph](https://cran.r-project.org/package=igraph) - Graph analysis for kinship matrices
- [future](https://cran.r-project.org/package=future) - Parallel processing framework
- [furrr](https://cran.r-project.org/package=furrr) - Future-based parallel mapping
- [R.utils](https://cran.r-project.org/package=R.utils) - Additional utility functions

### Running the Pipeline

1. **Sequential Execution**: Run steps 1-8 in order, as each step depends on outputs from previous steps.

2. **Path Configuration**: Each script has a standardized CONFIGURATION section at the top. Modify the placeholder paths for your environment:
   - Input data directories
   - Output directories  
   - Genotype data locations
   - Reference files

3. **No Additional Setup Required**: After configuring paths, scripts run autonomously without requiring further user input.


### Pipeline Steps

### Step 1: Data Processing
**Script**: [step1_process_olink_data.R](https://github.com/statgenetics/TransCisPredict/blob/main/scripts/step1_process_olink_data.R)

Processes raw OLINK proteomics data, performs quality control, and creates individual protein expression files.

### Step 2: Covariate Regression  
**Script**: [step2_covariate_regression.R](https://github.com/statgenetics/TransCisPredict/blob/main/scripts/step2_covariate_regression.R)

Removes covariate effects (age, sex, age×sex, BMI, genetic PCs 1-20) from protein expression to generate residuals for genetic analysis.

### Step 3: LD Block Selection
**Script**: [step3_LD_block_selection.R](https://github.com/statgenetics/TransCisPredict/blob/main/scripts/step3_LD_block_selection.R)

Identifies genomic regions with significant genetic signal using FDR-based selection to focus downstream analysis.

### Step 4: Cross Validation
**Script**: [step4_cross_validation.R](https://github.com/statgenetics/TransCisPredict/blob/main/scripts/step4_cross_validation.R)

Performs cross-validation to evaluate multiple statistical methods (BayesR, LASSO, Elastic Net, SuSiE) for protein prediction using both cis and trans genetic variants.

### Step 5: CV Evaluation
**Scripts**: 
- [step5a_evaluate_cv_performance.R](https://github.com/statgenetics/TransCisPredict/blob/main/scripts/step5a_evaluate_cv_performance.R) - Evaluates cross-validation performance
- [step5b_identify_best_method.R](https://github.com/statgenetics/TransCisPredict/blob/main/scripts/step5b_identify_best_method.R) - Selects optimal method per protein  
- [step5c_summarize_all_methods.R](https://github.com/statgenetics/TransCisPredict/blob/main/scripts/step5c_summarize_all_methods.R) - Comprehensive method summary (optional)

Evaluates cross-validation results and identifies the best-performing statistical method for each protein.

### Step 6: Whole Sample Analysis
**Script**: [step6_whole_sample_analysis.R](https://github.com/statgenetics/TransCisPredict/blob/main/scripts/step6_whole_sample_analysis.R)

Trains final prediction models using the complete dataset and optimal method identified for each protein. Generates final prediction weights for proteins meeting performance thresholds.

### Step 7: Population Prediction
**Scripts**:
- [step7a_predict_npx_population.R](https://github.com/statgenetics/TransCisPredict/blob/main/scripts/step7a_predict_npx_population.R) - Predicts protein levels for target population
- [step7b_combine_npx_files.R](https://github.com/statgenetics/TransCisPredict/blob/main/scripts/step7b_combine_npx_files.R) - Combines individual protein predictions into final matrix

Applies trained models to predict protein levels for UK Biobank individuals without measured proteomic data.

### Step 8: PWAS Analysis  
**Script**: [step8_pwas_analysis.R](https://github.com/statgenetics/TransCisPredict/blob/main/scripts/step8_pwas_analysis.R)

Performs proteome-wide association studies using generalized estimating equations to test associations between predicted proteins and complex traits.

### Utilities

**Directory**: [utilities/](https://github.com/statgenetics/TransCisPredict/tree/main/scripts/utilities)
**Scripts**:
- [pqtl_functions.R](https://github.com/statgenetics/TransCisPredict/blob/main/scripts/utilities/pqtl_functions.R) - Core utility functions for data processing and genetic prediction
- [pqtl_weights.R](https://github.com/statgenetics/TransCisPredict/blob/main/scripts/utilities/pqtl_weights.R) - Statistical method implementations (BayesR, LASSO, Elastic Net, SuSiE)
- [timing_function.R](https://github.com/statgenetics/TransCisPredict/blob/main/scripts/utilities/timing_function.R) - Runtime calculation utilities

Contains common functions and method implementations used throughout the pipeline.

## Key Features

- **Standardized Configuration**: All scripts use consistent CONFIGURATION sections for easy setup
- **Comprehensive Documentation**: Each parameter includes detailed input/output specifications
- **Flexible Methods**: Supports multiple statistical approaches with automatic method selection
- **Scalable Processing**: Designed for large-scale proteomics and genomics data
- **Quality Control**: Built-in validation and error handling throughout the pipeline

This pipeline provides a complete framework for protein expression imputation and proteome-wide association studies incorporating both cis and trans genetic variants for improved prediction performance.
