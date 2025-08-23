# ============================================================================
# Step 5c: Summarize All Cross-Validation Methods (Comprehensive Analysis)
# ============================================================================
# 
# Purpose:
# This script provides comprehensive performance metrics for ALL statistical 
# methods across all proteins. Unlike Step 5b which identifies the single best
# method per protein, this script computes detailed performance summaries for
# BayesR, LASSO, Elastic Net, and SuSiE methods to enable cross-method comparison.
#
# Input:
# 1. CV prediction accuracy files from Step 5a:
#    - Files: {protein}_cv_prediction_accuracy.csv
#    - Contains performance metrics for all methods across CV folds
#    - Format: [cv_num, method1_corr, method1_r2, ..., method2_corr, method2_r2, ...]
#    - Located in prediction_accuracy/ subdirectory from Step 5a output
#
# Output:
# 1. Comprehensive methods summary: all_methods_wide_by_cv_for_protein.csv
#    - Performance metrics for ALL methods across all proteins
#    - Format: [protein, method1_avg_r2, method1_sd_r2, method1_avg_corr, method2_avg_r2, ...]
#    - Contains mean R², standard deviation R², and mean correlation for each method
#    - One row per protein, columns for each method's performance metrics
#    - Enables comparison of method performance patterns across the proteome
#
# Usage Example:
# Rscript summarize_all_methods.R path_to_prediction_accuracy_files path_to_output_dir
#
# Command line arguments:
# Argument 1: input_dir - Directory containing *_cv_prediction_accuracy.csv files from Step 5a
# Argument 2: output_dir - Directory to write comprehensive methods summary file
#
# Prerequisites:
# - Must run AFTER Step 5a (evaluate_cv_performance.R)
# - Can be run independently of Step 5b (identify_best_method.R)
# - Useful for manuscript figures and method performance analysis
# ============================================================================

# Load required packages
suppressMessages({
    library(tidyverse)
})

## CONFIGURATION - MODIFY THESE PATHS FOR YOUR ENVIRONMENT

# Command line arguments
args = commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
    stop("Usage: Rscript summarize_all_methods.R <input_directory> <output_directory>")
}

input_dir = args[1]      # Directory with prediction accuracy files from Step 5a  
output_dir = args[2]     # Directory to write comprehensive methods summary

prediction_list <- list.files(input_dir, ".csv")

# Initialize output
output <- tibble()

for(file_name in prediction_list) {
    # Get protein name
    protein_name <- str_extract(file_name, "^[:graph:]+(?=_cv_)")

    data <- read_csv(paste0(input_dir, file_name), show_col_types = FALSE)

    # Identify unique methods (based on "_r2" columns)
    method_names <- colnames(data) |> 
        str_subset("_r2$") |> 
        str_remove("_r2$")

    # Collect summaries for all methods into one row
    method_summary <- map(method_names, function(method) {
        r2_col <- paste0(method, "_r2")
        corr_col <- paste0(method, "_corr")

        tibble(
            !!paste0(method, "_avg_r2")   := mean(data[[r2_col]], na.rm = TRUE),
            !!paste0(method, "_sd_r2")    := sd(data[[r2_col]], na.rm = TRUE),
            !!paste0(method, "_avg_corr") := mean(data[[corr_col]], na.rm = TRUE)
        )
    }) |> bind_cols()

    # Add protein column
    protein_summary <- tibble(protein = protein_name) |> bind_cols(method_summary)

    output <- bind_rows(output, protein_summary)

    print(paste0("Finished with protein: ", protein_name))
}

# Write output
write_csv(output, paste0(output_dir, "all_methods_wide_by_cv_for_protein.csv"))

