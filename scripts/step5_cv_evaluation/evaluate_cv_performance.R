# ============================================================================
# Step 5a: Cross-Validation Performance Evaluation 
# ============================================================================
# 
# Purpose:
# This script evaluates the performance of cross-validation results from Step 4.
# For each protein, it applies the trained weights from each CV fold to the 
# corresponding held-out test data and computes prediction accuracy metrics.
# This provides evaluation of how well each statistical method (BayesR, LASSO, 
# Elastic Net, SuSiE) performed during cross-validation.
#
# Input:
# 1. CV weights from Step 4:
#    - Files: {protein}_{region}_posterior_weights.csv
#    - Contains trained weights for all methods and CV folds
#    - Located in CV weights output directory from Step 4
#
# 2. Protein residuals from Step 2:  
#    - Files: {protein_name}_residual_npx.csv
#    - Covariate-adjusted NPX residuals for comparison with predictions
#    - Located in protein residuals directory from Step 2
#
# 3. Genotype data:
#    - PLINK format files (.bed/.bim/.fam) by LD blocks
#    - Same genotype data used in Steps 3-4
#    - Used to generate predictions for held-out test individuals
#
# Output:
# 1. Predicted NPX files: predicted_npx/{protein}_predicted_npx_from_covar_impute_missing_geno_cv{fold}.csv
#    - Contains predicted protein levels for each CV fold and statistical method
#    - Format: [eid, method1_weights, method2_weights, ...]
#    - One file per CV fold (cv1, cv2, cv3, cv4, cv5)
#
# 2. Accuracy metrics file: prediction_accuracy/{protein}_cv_prediction_accuracy.csv
#    - Performance metrics: correlation, R², adjusted R², RMSE, MAE, p-value
#    - Format: [cv_num, method1_corr, method1_r2, ..., method2_corr, method2_r2, ...]
#    - One row per CV fold, columns for each method's performance metrics
#
# Usage Example (for a1bg protein):
# Rscript evaluate_cv_performance.R a1bg path_to_step4_cv_weights path_to_step5_output
#
# Command line arguments:
# Argument 1: protein_name - Name of protein to evaluate (e.g., "a1bg", "APOE")
# Argument 2: input_dir - Directory containing CV weights from Step 4 
# Argument 3: output_dir - Directory to write evaluation results
#
# Prerequisites:
# - Must run AFTER Step 4 (cross-validation analysis)
# - Must run BEFORE Step 5b (identify_best_method.R) and 5c (summarize_all_methods.R)
# ============================================================================

## CONFIGURATION - MODIFY THESE PATHS FOR YOUR ENVIRONMENT

# Command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
    stop("Usage: Rscript evaluate_cv_performance.R <protein_name> <weights_directory> <output_directory>")
}

protein_name <- args[1]           # Protein to evaluate (e.g., "APOE")
input_dir <- args[2]              # Directory with CV weights from Step 4
output_dir <- args[3]             # Directory to write results

# Input paths - modify for your environment  
protein_residuals_dir <- "path_to_save_protein_residuals"     # Step 2 output
genotype_base_path <- "path_to_genotype_data_by_ld_blocks"    # Genotype data

# Analysis parameters
num_folds <- 5                    # Number of CV folds (must match Step 4)

# Create output directories if they don't exist
dir.create(file.path(output_dir, "predicted_npx"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_dir, "prediction_accuracy"), showWarnings = FALSE, recursive = TRUE)

# Start time to later measure total run time
time_1 <- Sys.time()

# Load packages
suppressMessages({
    library(tidyverse)
    library(data.table)
    library(plink2R)
    library(rsample)
    source("../utilities/pqtl_functions.R")
    source("../utilities/timing_function.R")
})

## PART 2: LOAD PROTEIN RESIDUALS
cat("Loading protein residuals for", protein_name, "...\n")

# Load protein residuals from Step 2
protein_residual_file <- file.path(protein_residuals_dir, paste0(protein_name, "_residual_npx.csv"))

if (!file.exists(protein_residual_file)) {
    stop("Protein residual file not found: ", protein_residual_file)
}

# Load protein residuals
protein_levels <- fread(protein_residual_file)

# Read in indvidual data
protein_levels <- protein_levels |>
    rename_with(~"participant_id", .cols = 1) |>
    rename_with(~"npx_res", .cols = 2) |>
    arrange(participant_id) |>
    drop_na() 

# Find CV weight files for this protein
post_weights_list <- list.files(input_dir, pattern = paste0(protein_name, ".*posterior_weights\\.csv$"))

if (length(post_weights_list) == 0) {
    stop("No weight files found for protein ", protein_name, " in directory: ", input_dir)
}

cat("Found", length(post_weights_list), "weight files for", protein_name, "\n")

# Initialize output list
results <- list()

for (cv_num in 1:num_folds){
    results[[cv_num]] <- 0
}

## PART 3: REGION CONTRIBUTION
# Loop over the directory of posterior weights
for (file in post_weights_list){

    # Start region
    print(paste0("Starting region ", str_extract(file, paste0("(?<=", protein_name, "_)[:graph:]+(?=_post)"))))
    chr_char <- str_extract(file, paste0("(?<=", protein_name,"_)[:alnum:]+(?=_)"))
    
    # Load in and prepare weights
    weights <- fread(paste0(input_dir, file)) |>
        drop_na()
    
    weights_matrix <- as.matrix(weights[, -1])
    rownames(weights_matrix) <- weights$variant_id
    
    if (any(str_detect(colnames(weights_matrix), "lm_weights"))) {
        next
    }
    
    # Load genotype data for this LD block
    region_id <- str_extract(file, paste0("(?<=", protein_name,"_)[:graph:]+(?=_post)"))
    genotype_file <- file.path(genotype_base_path, chr_char, region_id)
    
    if (!file.exists(paste0(genotype_file, ".bed"))) {
        cat("Warning: Genotype file not found:", paste0(genotype_file, ".bed"), "- skipping region\n")
        next
    }
    
    plink2r_data <- read_plink(genotype_file)
    
    # Extract genotype matrix
    gen_matrix <- plink2r_data$bed
    rownames(gen_matrix) <- str_extract(rownames(gen_matrix), "(?<=:)(\\d+)")
    
    # QC the rows
    common_ids <- intersect(protein_levels$participant_id, rownames(gen_matrix))
    gen_matrix <- gen_matrix[common_ids, , drop = FALSE]
    
    # QC the columns
    rsid_overlap = intersect(rownames(weights_matrix), colnames(gen_matrix))
    if (length(rsid_overlap) == 0) stop("Misalignment between weights and genotype")
    
    weights_matrix <- weights_matrix[rsid_overlap,, drop=FALSE]
    gen_matrix <- gen_matrix[,rsid_overlap, drop=FALSE]
    
    gen_matrix_imputed <- apply(gen_matrix, 2, function(col) {
    if (any(is.na(col))) {
        col[is.na(col)] <- stat_mode(col[!is.na(col)])
    }
        return(col)
    })
    
    # Perform CV. Set seed for reproducibility
    set.seed(1)
    df_cv <- vfold_cv(gen_matrix_imputed, v = num_folds) 
    
    for (cv_num in 1:num_folds){
        mat_cv <- assessment(df_cv[[1]][[cv_num]])
        
        # select only the weights for this fold
        mat_weights <- colnames(weights_matrix) %>%
            str_detect(paste0("_", cv_num)) %>%
            weights_matrix[,.]
        
        results[[cv_num]] <- results[[cv_num]] + (mat_cv %*% mat_weights)
        colnames(results[[cv_num]]) <- colnames(mat_weights)
        
    }
}


## PART 4: PREDICTION ACCURACY
# Initialize output of summary
prediction_summary <- tibble()

# Loop over cv folds
for (cv_num in 1:num_folds){
    # Subset the protein npx_res
    observed_protein <- protein_levels |>
        filter(participant_id %in% rownames(results[[cv_num]]))
    
    fold_summary <- c("cv_num" = cv_num)
    
    for (j in 1:ncol(results[[cv_num]])){
        # Get name and predicted expression for method
        method <- colnames(results[[cv_num]]) |>
            nth(j) |>
            str_extract("^[:graph:]+(?=_weights)")

        predicted_protein <- results[[cv_num]][,j]
        
        if (all(predicted_protein == 0)) {
            fold_corr <- fold_r2 <- fold_adj_r2 <- 0
            fold_rmse <- fold_mae <- fold_pval <- NA
        } else {
        
            # Fit lm
            fit <- lm(observed_protein$npx_res ~ predicted_protein)
        
            # Get performance output
            fold_corr <- cor(observed_protein$npx_res, predicted_protein)
            fold_r2 <- summary(fit)$r.squared
            fold_adj_r2 <-summary(fit)$adj.r.squared
            fold_rmse <- (observed_protein$npx_res - predicted_protein)^2 |>
                mean() |>
                sqrt()
            fold_mae <- abs(observed_protein$npx_res - predicted_protein) |>
                mean()     
            fold_pval <- summary(fit)$coefficients[2,4]
        }
        
        # Update output
        method_summary <- c(fold_corr, fold_r2, fold_adj_r2, fold_rmse, fold_mae, fold_pval)
        names(method_summary) <- c(paste0(method, c("_corr", "_r2", "_adj_r2", 
                                                    "_rmse", "_mae", "_pval")))
        fold_summary <- c(fold_summary, method_summary)
    }

    
    # Add cv to tibble
    prediction_summary <- rbind(prediction_summary, fold_summary)
    colnames(prediction_summary) <- c(names(fold_summary))
}

## PART 5: OUTPUT
# Write outputs of prediction
for (cv_num in 1:num_folds){
    # Write output
    results[[cv_num]] |>
        as_tibble() |>
        mutate(eid = rownames(results[[cv_num]])) |>
        write_csv(file = paste0(output_dir, "predicted_npx/",
                            protein_name, "_predicted_npx_from_covar_impute_missing_geno_cv", cv_num, ".csv"))
}

# Write accuracy metrics
prediction_summary |>
    write_csv(file = paste0(output_dir, "prediction_accuracy/", protein_name, "_cv_prediction_accuracy.csv"))


## PART 6: COMPLETION SUMMARY
time_2 <- Sys.time()
total_duration <- hms_span(time_1, time_2)

cat("\n============================================================================\n")
cat("Cross-Validation Evaluation Completed!\n")
cat("============================================================================\n")
cat("Protein:", protein_name, "\n")
cat("Weight files processed:", length(post_weights_list), "\n")
cat("CV folds evaluated:", num_folds, "\n")
cat("Output directory:", output_dir, "\n")
cat("Total runtime:", total_duration, "\n")
cat("============================================================================\n")
cat("Next steps:\n")
cat("1. Run identify_best_method.R to find best performing method\n") 
cat("2. OR run summarize_all_methods.R to get comprehensive performance summary\n")
cat("============================================================================\n")
