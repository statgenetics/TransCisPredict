# ============================================================================
# Step 5b: Identify Best Cross-Validation Method Per Protein
# ============================================================================
# Identifies the best-performing statistical method for each protein based on 
# cross-validation results from Step 5a. Compares performance across all methods 
# and selects the method with highest average R² across CV folds.
#
# Prerequisites: Must run AFTER Step 5a

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

# OUTPUT: Best methods summary
# - File: best_method_by_cv_for_protein.csv
# - Contains best method identification for each protein
# - Format: [protein, best_method, best_avg_r2, avg_cv_corr, best_r2_sd]
# - best_method: Method name with highest average R² (e.g., "bayesr_weights")
# - best_avg_r2: Average R² across CV folds for best method
# - avg_cv_corr: Average correlation across CV folds for best method  
# - best_r2_sd: Standard deviation of R² across CV folds for best method
output_dir <- "path_to_step5b_output/"

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
    
    data = read_csv(paste0(input_dir, file_name), show_col_types=FALSE)
    # Load predictions
    df_pred <- read_csv(paste0(input_dir, file_name), show_col_types=FALSE) |>
        select(-ends_with("adj_r2")) |>
        select(ends_with("r2"))
    
    # Identify best method
    best_method <- df_pred |>
        colMeans() |>
        which.max() |>
        names()
    
    # Best average r^2
    best_r2 <- df_pred |>
        colMeans() |>
        max()
    
    avg_cv_corr <- data |>
        pull(paste0(str_remove(best_method, "_r2$"), "_corr")) |>
        mean()

    best_r2_sd <- df_pred |>
        pull(best_method) |>
        sd()
    
    # Update name
    best_method <- str_replace(best_method, "r2", "weights")
    
    # Add protein to list
    protein_summary <- tibble(protein = protein_name, 
                              best_method = best_method,
                              avg_r2 = best_r2,
                              sd_r2 = best_r2_sd,
                              avg_corr = avg_cv_corr
)
    
    output <- rbind(output, protein_summary)
    
    print(paste0("Finished with protein: ", protein_name))
}

# Write output
write_csv(output, paste0(output_dir, "best_method_by_cv_for_protein.csv"))
