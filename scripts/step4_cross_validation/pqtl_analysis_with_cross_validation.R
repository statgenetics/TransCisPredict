# ============================================================================
# Step 4: Cross-Validation Analysis for Protein Prediction Methods
# ============================================================================
# 
# Purpose:
# This script performs 5-fold cross-validation to evaluate different statistical 
# methods for predicting protein levels from genetic data. For each protein and
# genomic region, it trains multiple methods (BayesR, LASSO, Elastic Net, SuSiE) 
# and generates prediction weights for performance comparison in Step 5.
# Optionally generates Step 5 job scripts for batch processing systems.
#
# Input:
# 1. Protein residuals from Step 2:
#    - Files: {protein_name}_residual_npx.csv
#    - Format: [eid, npx_residuals]
#    - Covariate-adjusted protein phenotypes
#
# 2. FDR region summaries from Step 3:
#    - Files: {protein_name}_region_stats_summary.csv
#    - Used to filter regions with genetic signal (signal_FDR0.2 = 1)
#
# 3. Genotype data by LD blocks:
#    - PLINK format files (.bed/.bim/.fam)
#    - Directory structure: chromosome/LD_block files
#    - Same data as used in Step 3
#
# Output:
# CV prediction weights: {protein}_{region}_posterior_weights.csv
# - Columns: variant_id, {method1}_1, {method1}_2, ..., {method4}_5
# - 4 methods × 5 CV folds = 20 weight columns per variant
# - One file per genomic region with significant signal
# - Used in Step 5 for performance evaluation
#
# Methods Evaluated:
# - BayesR: Bayesian variable selection regression
# - LASSO: L1-regularized regression
# - Elastic Net: Combined L1/L2 regularization  
# - SuSiE: Sum of Single Effects fine-mapping
# (Linear regression used for single-variant regions)
#
# Processing:
# - 5-fold cross-validation with consistent random seeds
# - Only analyze regions with significant genetic signal from Step 3
# - Train each method on 4/5 of data, generate weights for prediction
# - Save weights for all methods and CV folds
# ============================================================================

# Load required packages and functions
suppressMessages({
    library(tidyverse)
    library(rsample)      # For cross-validation splitting
    library(plink2R)      # For genotype data loading
    library(pecotmr)      # For statistical methods (BayesR, SuSiE, etc.)
    library(data.table)   # For efficient data handling
    library(future)       # For parallel processing
    library(furrr)        # For parallel map functions
    library(janitor)      # For column name cleaning
    source("../utilities/pqtl_weights.R")    # Statistical methods implementation
    source("../utilities/pqtl_functions.R")  # Utility functions
    source("../utilities/timing_function.R") # Timing utilities
})

# ============================================================================
# CONFIGURATION - MODIFY THESE PATHS FOR YOUR ENVIRONMENT  
# ============================================================================

# Command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
    stop("Usage: Rscript pqtl_analysis_with_cross_validation.R <protein_name> <chromosome> <output_directory> <fdr_threshold> [generate_step5_scripts]")
}

protein_name <- args[1]           # Protein to analyze
chr_num <- args[2]               # Chromosome (e.g., "01", "02")
output_dir <- args[3]            # Output directory for weights
min_q_threshold <- args[4]       # FDR threshold (0.2 for main results)
generate_step5_scripts <- if(length(args) >= 5) as.logical(args[5]) else FALSE  # Optional: generate Step 5 job scripts

if (is.na(min_q_threshold)) min_q_threshold <- 0.2

# Analysis parameters
num_folds <- 5                   # Number of CV folds
cv_seed <- 1                     # Random seed for reproducibility

# Input paths - modify for your environment
protein_residuals_dir <- "path_to_save_protein_residuals"     # Step 2 output
fdr_summaries_dir <- "path_to_fdr_region_summaries"          # Step 3 output  
genotype_base_path <- "path_to_genotype_data_by_ld_blocks"   # Same as Step 3

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}

# Start timing
time_1 <- Sys.time()

cat("============================================================================\n")
cat("Cross-Validation Analysis\n") 
cat("============================================================================\n")
cat("Protein:", protein_name, "\n")
cat("Chromosome:", chr_num, "\n")
cat("FDR threshold:", min_q_threshold, "\n")
cat("Number of CV folds:", num_folds, "\n")
cat("Output directory:", output_dir, "\n")
cat("Generate Step 5 scripts:", generate_step5_scripts, "\n")
cat("============================================================================\n")

# ============================================================================
# DATA LOADING
# ============================================================================

cat("\nLoading input data...\n")

# Load protein residuals from Step 2
protein_residual_file <- file.path(protein_residuals_dir, paste0(protein_name, "_residual_npx.csv"))

if (!file.exists(protein_residual_file)) {
    stop("Protein residual file not found: ", protein_residual_file)
}

protein_levels <- fread(protein_residual_file) %>%
    rename_with(~"participant_id", .cols = 1) %>%
    rename_with(~"npx", .cols = 2) %>%
    arrange(participant_id) %>%
    drop_na()

cat("Loaded protein residuals for", nrow(protein_levels), "individuals\n")

# Load FDR summary from Step 3 to identify regions with signal
fdr_summary_file <- file.path(fdr_summaries_dir, paste0(protein_name, "_region_stats_summary.csv"))

if (!file.exists(fdr_summary_file)) {
    stop("FDR summary file not found: ", fdr_summary_file)
}

fdr_summary <- fread(fdr_summary_file)
cat("Loaded FDR summary for", nrow(fdr_summary), "regions\n")

# Get list of regions for this chromosome
chr_genotype_dir <- file.path(genotype_base_path, chr_num)
if (!dir.exists(chr_genotype_dir)) {
    stop("Chromosome genotype directory not found: ", chr_genotype_dir)
}

region_list <- list.files(
    path = chr_genotype_dir,
    pattern = "\\.bed$",  # Look for .bed files (PLINK binary format)
    full.names = TRUE
)

if (length(region_list) == 0) {
    stop("No .bed files found in chromosome directory: ", chr_genotype_dir)
}

cat("Found", length(region_list), "LD regions in chromosome", chr_num, "\n")

# ============================================================================
# CROSS-VALIDATION ANALYSIS BY REGION
# ============================================================================

regions_processed <- 0
regions_skipped <- 0

# Process each genomic region
for (region_file in region_list) {
    
    # Extract region identifier from filename
    region_basename <- tools::file_path_sans_ext(basename(region_file))
    region_id <- str_extract(region_basename, "^[^_]+(?=_un|$)")
    
    cat("\n=========================", "Region", region_id, "==========================\n")
    cat("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
    
    # Check if region has significant genetic signal from Step 3
    fdr_column_name <- paste0("signal_FDR", min_q_threshold)
    region_row <- fdr_summary[fdr_summary$region_name == region_id, ]
    
    if (nrow(region_row) == 0) {
        cat("WARNING: Region", region_id, "not found in FDR summary - skipping\n")
        regions_skipped <- regions_skipped + 1
        next
    }
    
    region_has_signal <- region_row[[fdr_column_name]]
    min_pval <- region_row$min_pval
    
    cat("Minimum p-value:", min_pval, "\n")
    
    # Skip regions without significant signal (flat regions)
    if (is.na(region_has_signal) || region_has_signal == 0) {
        cat("Region classification: FLAT (no significant signal) - skipping\n")
        regions_skipped <- regions_skipped + 1
        next
    }
    
    cat("Region classification: SIGNAL detected - proceeding with analysis\n")
    
    tryCatch({
        
        # Load genotype data for this region
        import_start <- Sys.time()
        plink_path <- tools::file_path_sans_ext(region_file)
        plink2r_data <- read_plink(plink_path)
        gen_matrix <- plink2r_data$bed
        
        # Extract participant IDs from rownames (format: "famid:iid")
        rownames(gen_matrix) <- str_extract(rownames(gen_matrix), "(?<=:)(\\d+)")
        rm(plink2r_data)  # Free memory
        import_end <- Sys.time()
        
        # Find common individuals between protein and genotype data
        merge_start <- Sys.time()
        common_ids <- intersect(protein_levels$participant_id, as.numeric(rownames(gen_matrix)))
        
        if (length(common_ids) < 100) {
            cat("Warning: Only", length(common_ids), "overlapping individuals - skipping region\n")
            regions_skipped <- regions_skipped + 1
            next
        }
        
        # Filter both datasets to common individuals
        protein_levels_matched <- protein_levels %>% filter(participant_id %in% common_ids)
        gen_matrix <- gen_matrix[as.character(common_ids), , drop = FALSE]
        merge_end <- Sys.time()
        
        # Report timing
        import_duration <- import_end - import_start
        merge_duration <- merge_end - merge_start
        cat("Step 1. Import genotype data:", format(import_duration, digits = 3), "\n")
        cat("Step 2. Merge phenotype/genotype data:", format(merge_duration, digits = 3), "\n")
        
        # Start main analysis
        pqtl_analysis_start <- Sys.time()
        
        # Combine protein and genotype data
        df_analysis <- cbind(protein_levels_matched, gen_matrix) %>%
            clean_names()  # Standardize column names, handle duplicate rsIDs
        
        cat("Analysis dataset:", nrow(df_analysis), "individuals ×", ncol(gen_matrix), "variants\n")
        
        # Configure statistical methods based on number of variants
        if (ncol(gen_matrix) == 1) {
            # Use simple linear regression for single-variant regions
            weight_methods <- list(lm_weights = list())
            cat("Single variant region - using linear regression\n")
        } else {
            # Use full set of methods for multi-variant regions  
            weight_methods <- list(
                bayes_r_weights = list(nit = 1000, nburn = 200, nthin = 5),  # MCMC parameters
                susie_weights = list(),
                lasso_weights = list(),
                enet_weights = list()
            )
            cat("Multi-variant region - using 4 methods:", paste(names(weight_methods), collapse = ", "), "\n")
        }
        
        # Perform 5-fold cross-validation
        set.seed(cv_seed)  # Ensure reproducibility
        df_cv <- vfold_cv(df_analysis, v = num_folds)
        cat("Created", num_folds, "-fold cross-validation splits\n")
        
        # Initialize storage for all CV weights
        post_weights <- list()
        
        # Train models on each CV fold
        for (cv_num in 1:num_folds) {
            cat("  Training fold", cv_num, "/", num_folds, "...")
            
            # Extract training data for this fold
            df_train <- analysis(df_cv$splits[[cv_num]])
            
            # Prepare phenotype and genotype matrices
            y_train <- df_train$npx
            X_train <- df_train %>%
                select(-participant_id, -npx) %>%
                as.matrix()
            
            # Run statistical methods to generate weights
            fold_weights <- pqtl_weights(
                y = y_train,
                X = X_train, 
                weight_methods = weight_methods,
                num_threads = 1
            )
            
            post_weights[[cv_num]] <- fold_weights
            cat("completed\n")
        }
        
        # Format output weights matrix
        num_methods <- length(weight_methods)
        methods <- names(weight_methods)
        
        # Create column names: method1_fold1, method1_fold2, ..., method4_fold5
        col_names <- paste0(rep(methods, each = num_folds), "_", rep(1:num_folds, times = num_methods))
        snp_names <- colnames(df_train)[-c(1, 2)]  # Exclude participant_id and npx
        
        # Initialize output matrix
        output <- matrix(NA, nrow = length(snp_names), ncol = length(col_names))
        rownames(output) <- snp_names
        colnames(output) <- col_names
        output <- as.data.frame(output)
        
        # Fill output matrix with weights from each fold and method
        for (i in 1:num_folds) {
            fold_weights <- post_weights[[i]]
            
            for (method in methods) {
                if (method %in% names(fold_weights)) {
                    weights <- fold_weights[[method]]
                    
                    if (!is.null(weights) && nrow(weights) > 0) {
                        weight_snp_ids <- rownames(weights)
                        output_col_name <- paste0(method, "_", i)
                        
                        # Match SNP IDs and fill corresponding weights
                        matched_indices <- match(weight_snp_ids, rownames(output))
                        valid_indices <- !is.na(matched_indices)
                        
                        if (sum(valid_indices) > 0) {
                            output[matched_indices[valid_indices], output_col_name] <- weights[valid_indices, 1]
                        }
                    }
                }
            }
        }
        
        # Add variant ID column and save results
        output_final <- output %>%
            as_tibble() %>%
            mutate(variant_id = rownames(output)) %>%
            select(variant_id, everything())
        
        # Save CV weights to output file
        output_file <- file.path(output_dir, paste0(protein_name, "_", region_id, "_posterior_weights.csv"))
        write_csv(output_final, output_file)
        
        pqtl_analysis_end <- Sys.time()
        pqtl_analysis_duration <- pqtl_analysis_end - pqtl_analysis_start
        
        cat("Step 3. Cross-validation analysis for", length(common_ids), "individuals took:", 
            format(pqtl_analysis_duration, digits = 3), "\n")
        cat("Saved weights to:", basename(output_file), "\n")
        
        # Clean up memory
        rm(gen_matrix, post_weights, protein_levels_matched, df_analysis, output, output_final)
        
        regions_processed <- regions_processed + 1
        
    }, error = function(e) {
        cat("Error processing region", region_id, ":", e$message, "\n")
        regions_skipped <- regions_skipped + 1
    })
}

# ============================================================================
# OPTIONAL: GENERATE STEP 5 JOB SCRIPTS
# ============================================================================

if (generate_step5_scripts && regions_processed > 0) {
    cat("\n============================================================================\n")
    cat("Generating Step 5 Job Scripts\n")
    cat("============================================================================\n")
    
    # Configuration for Step 5 script generation
    step5_script_dir <- "path_to_step5_scripts_directory"
    step5_output_dir <- "path_to_step5_output_directory" 
    step5_log_dir <- "path_to_step5_logs_directory"
    step4_weights_dir <- output_dir  # Use current output directory as input to Step 5
    
    # Create directories if they don't exist
    for (dir_path in c(step5_script_dir, step5_output_dir, step5_log_dir)) {
        if (!dir.exists(dir_path)) {
            dir.create(dir_path, recursive = TRUE)
        }
    }
    
    # Create Step 5 subdirectories
    dir.create(file.path(step5_output_dir, "predicted_npx"), showWarnings = FALSE, recursive = TRUE)
    dir.create(file.path(step5_output_dir, "prediction_accuracy"), showWarnings = FALSE, recursive = TRUE)
    
    # Generate SLURM job script for Step 5 analysis
    job_script_path <- file.path(step5_script_dir, paste0(protein_name, "_step5_prediction.sbatch"))
    
    job_script_content <- paste0(
        "#!/bin/bash\n",
        "#SBATCH --job-name=predict_npx_", protein_name, "\n",
        "#SBATCH --mem=15G\n", 
        "#SBATCH --time=24:00:00\n",
        "#SBATCH --output=", step5_log_dir, "/", protein_name, "_npx_prediction_%j.out\n",
        "#SBATCH --error=", step5_log_dir, "/", protein_name, "_npx_prediction_%j.err\n",
        "\n",
        "mkdir -p ", step5_output_dir, "\n",
        "\n",
        "# Run Step 5 prediction analysis using Step 4 cross-validation weights\n",
        "Rscript path_to_step5_script/predict_npx_from_cv.R ", protein_name, " ", step4_weights_dir, "/ ", step5_output_dir, 
        " &> ", step5_log_dir, "/", protein_name, "_npx_prediction.log\n"
    )
    
    # Write job script to file
    writeLines(job_script_content, job_script_path)
    
    cat("Generated Step 5 job script:", basename(job_script_path), "\n")
    cat("Step 4 weights directory (input to Step 5):", step4_weights_dir, "\n")
    cat("Step 5 output directory:", step5_output_dir, "\n")
    cat("Step 5 logs directory:", step5_log_dir, "\n")
}

# ============================================================================
# SUMMARY AND COMPLETION
# ============================================================================

time_2 <- Sys.time()
total_duration <- hms_span(time_1, time_2)

cat("\n============================================================================\n")
cat("Cross-Validation Analysis Completed!\n")
cat("============================================================================\n")
cat("Protein:", protein_name, "\n")
cat("Chromosome:", chr_num, "\n")
cat("Regions processed:", regions_processed, "\n")
cat("Regions skipped:", regions_skipped, "(no signal or errors)\n")
cat("Total runtime:", total_duration, "\n")
cat("Output directory:", output_dir, "\n")
cat("============================================================================\n")

# ============================================================================
# EXPECTED OUTPUT FORMAT
# ============================================================================
#
# Output files: {protein}_{region}_posterior_weights.csv
# 
# Format:
# variant_id, bayes_r_weights_1, bayes_r_weights_2, ..., bayes_r_weights_5,
#             susie_weights_1, susie_weights_2, ..., susie_weights_5,
#             lasso_weights_1, lasso_weights_2, ..., lasso_weights_5,
#             enet_weights_1, enet_weights_2, ..., enet_weights_5
#
# - Each row represents a genetic variant
# - Each method has 5 columns (one per CV fold)
# - Total: 20 columns (4 methods × 5 folds) + variant_id
# - NA values indicate method did not select that variant
# - Non-NA values are prediction weights
#
# Usage in Step 5:
# - These weights are applied to test set genotypes
# - Performance metrics calculated for each method
# - Best performing method identified per protein
#
# Optional Step 5 Script Generation:
# - If generate_step5_scripts=TRUE, creates SLURM job scripts for Step 5
# - Generated scripts use Step 4 weights as input for protein prediction
# - Parameterized paths must be configured for your batch system environment
# ============================================================================