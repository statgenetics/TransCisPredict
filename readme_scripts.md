# TransCisPredict Analysis Pipeline

This directory contains the complete analysis pipeline for the TransCisPredict framework - a comprehensive approach to protein expression weights estimation that incorporates both _cis_- and _trans_- variants for performing proteome-wide association studies (PWAS) with predicted protein expression levels. The pipeline consists of eight sequential steps, each implemented as a standalone R script with standardized configuration sections.

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
**Script**: [step1_process_olink_data.R](https://github.com/statgenetics/TransCisPredict/blob/main/step1_process_olink_data.R)

Process raw protein expression levels, perform quality control, and create individual protein expression files.

### Step 2: Covariate Regression
**Script**: [step2_covariate_regression.R](https://github.com/statgenetics/TransCisPredict/blob/main/step2_covariate_regression.R)

Remove covariate effects, e.g., age, sex, age×sex, body mass index (BMI), genetic ancestry principal components (PCs) 1-20, from protein expression to generate normalized protein expression (NPX) residuals for each protein.

### Step 3: LD Block Selection
**Script**: [step3_LD_block_selection.R](https://github.com/statgenetics/TransCisPredict/blob/main/step3_LD_block_selection.R)

Identify genomic regions with significant genetic signal using false discovery rate (FDR)-based selection to include in weight estimation.

### Step 4: Cross Validation
**Script**: [step4_cross_validation.R](https://github.com/statgenetics/TransCisPredict/blob/main/step4_cross_validation.R)

Perform cross-validation to evaluate multiple statistical methods (BayesR, SuSiE, LASSO, and Elastic Net) for protein prediction using both _cis_- and _trans_- variants.

### Step 5: CV Evaluation
**Scripts**:
- [step5a_evaluate_cv_performance.R](https://github.com/statgenetics/TransCisPredict/blob/main/step5a_evaluate_cv_performance.R) - Evaluate cross-validation performance for four methods per protein
- [step5b_identify_best_method.R](https://github.com/statgenetics/TransCisPredict/blob/main/step5b_identify_best_method.R) - Select optimal method for each protein
- [step5c_summarize_all_methods.R](https://github.com/statgenetics/TransCisPredict/blob/main/step5c_summarize_all_methods.R) - Generate a comprehensive results summary for all proteins (optional), e.g., $CV-r$ and $CV-R^2$

Evaluate cross-validation results and identifies the "optimal" method for each protein.

### Step 6: Whole Sample Analysis
**Script**: [step6_whole_sample_analysis.R](https://github.com/statgenetics/TransCisPredict/blob/main/step6_whole_sample_analysis.R)

Estimate the weights using the complete reference sample applying the "optimal" method identified for each protein.

### Step 7: Prediction of NPX Residual Levels in Target Sample
**Scripts**:
- [step7a_predict_npx_population.R](https://github.com/statgenetics/TransCisPredict/blob/main/step7a_predict_npx_population.R) - Predict NPX residuals for each protein for every individual in the target sample
- [step7b_combine_npx_files.R](https://github.com/statgenetics/TransCisPredict/blob/main/step7b_combine_npx_files.R) - Combine all proteins predicted NPX residual levels into a single file

Apply weights to predict NPX residual levels in the target sample.

### Step 8: PWAS Analysis
**Script**: [step8_pwas_analysis.R](https://github.com/statgenetics/TransCisPredict/blob/main/step8_pwas_analysis.R)

Perform proteome-wide association studies using generalized estimating equations (GEE) to test associations between predicted NPX residual levels and complex traits.

### Utilities

**Directory**: [utilities/](https://github.com/statgenetics/TransCisPredict/tree/main/utilities)
**Scripts**:
- [pqtl_functions.R](https://github.com/statgenetics/TransCisPredict/blob/main/utilities/pqtl_functions.R) - Core utility functions for data processing and genetic prediction
- [pqtl_weights.R](https://github.com/statgenetics/TransCisPredict/blob/main/utilities/pqtl_weights.R) - Statistical method implementations (BayesR, SuSiE, LASSO, and Elastic Net)
- [timing_function.R](https://github.com/statgenetics/TransCisPredict/blob/main/utilities/timing_function.R) - Runtime calculation utilities

Contain common functions and method implementations used throughout the pipeline.
