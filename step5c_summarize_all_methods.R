# ============================================================================
# Step 5c: Summarize All Cross-Validation Methods (Comprehensive Analysis)
# ============================================================================
# Provides comprehensive performance metrics for ALL statistical methods across 
# all proteins. Unlike Step 5b which identifies the single best method per protein, 
# this script computes detailed performance summaries for BayesR, LASSO, Elastic 
# Net, and SuSiE methods to enable cross-method comparison.
# 
# NOTE: This script is OPTIONAL. Step 5b already identifies the best method for
# each protein. Use this script only if you need detailed performance comparisons
# across all methods for analysis or manuscript figures.

# Load packages
suppressMessages({
    library(tidyverse)
})

# ============================================================================
# CONFIGURATION - MODIFY THESE PATHS FOR YOUR ENVIRONMENT
# ============================================================================

# INPUT: CV prediction accuracy files from Step 5a
# - Files: {protein}_cv_prediction_accuracy.csv
# - Contains performance metrics for all methods across CV folds
# - Format: [cv_num, method1_corr, method1_r2, ..., method2_corr, method2_r2, ...]
# - Located in prediction_accuracy/ subdirectory from Step 5a output
input_dir <- "path_to_step5a_prediction_accuracy_files/"

# OUTPUT: Comprehensive methods summary
# - File: all_methods_wide_by_cv_for_protein.csv
# - Performance metrics for ALL methods across all proteins
# - Format: [protein, method1_avg_r2, method1_sd_r2, method1_avg_corr, method2_avg_r2, ...]
# - Contains mean R², standard deviation R², and mean correlation for each method
# - One row per protein, columns for each method's performance metrics
# - Enables comparison of method performance patterns across the proteome
output_dir <- "path_to_step5c_output/"

# END CONFIGURATION
# ============================================================================

# ============================================================================
# SCRIPT EXECUTION
# ============================================================================

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

