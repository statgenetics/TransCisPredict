# ============================================================================
# Step 7b: Combine Population Protein Predictions 
# ============================================================================
# Combines individual protein prediction files from Step 7a into a single 
# comprehensive dataset containing predicted protein levels for all proteins 
# across the entire UK Biobank White British cohort. Creates the final 
# proteome-wide prediction matrix for PWAS analysis.
#
# Prerequisites: Must run AFTER Step 7a

# Load packages
suppressMessages({
    library(tidyverse)
    library(data.table)
})

# ============================================================================
# CONFIGURATION - MODIFY THESE PATHS FOR YOUR ENVIRONMENT
# ============================================================================

# INPUT: Individual protein prediction files from Step 7a
# - Files: {protein}_predicted_npx_all_white_british.csv
# - Each contains: [participant_id, predicted_npx] for one protein
input_dir <- "path_to_step7a_individual_predictions"

# OUTPUT: Combined proteome predictions
# - File: combined_predicted_proteome_all_white_british.csv
# - Contains predicted NPX values for ALL proteins and ALL individuals
# - Format: [IID, protein1, protein2, protein3, ..., proteinN]
# - One row per individual, one column per protein
# - Used as input for Step 8 PWAS analysis
output_file <- "path_to_combined_proteome_predictions.csv"   # Full path to combined output file

# ============================================================================
# SCRIPT EXECUTION
# ============================================================================

cat("Starting protein prediction file combination...\n")
cat("Input directory:", input_dir, "\n")
cat("Output file:", output_file, "\n")

# Get list of prediction files
file_list <- list.files(input_dir, pattern = "_predicted_npx_.*\\.csv$", full.names = TRUE)

if (length(file_list) == 0) {
    stop("No prediction files found in directory: ", input_dir)
}

cat("Found", length(file_list), "protein prediction files\n")

# Initialize output with participant IDs from first file
cat("Initializing combined dataset...\n")
output = fread(file_list[1]) |>
    select(IID)

expected_rows <- nrow(output)
cat("Expected number of individuals per file:", expected_rows, "\n")

# Combine all protein predictions
for (file in file_list){
    # Extract protein name from filename
    protein_name = str_extract(basename(file), "^[:graph:]+(?=_predicted)")
    
    # Read and process protein predictions
    df = fread(file) |>
        select(-any_of("FID")) |>  # Remove FID if present
        rename(!!sym(protein_name) := npx)
    
    # Quality check: verify row count matches expected
    if (nrow(df) != expected_rows) {
        stop("Error: ", basename(file), " has ", nrow(df), " rows, expected ", expected_rows)
    }
    
    # Merge with main output
    output = left_join(output, df, by = join_by(IID))
    
    # Progress update
    cat("Processed protein", which(file_list == file), "of", length(file_list), ":", protein_name, "\n")
}

# Write combined output
cat("Writing combined proteome predictions to:", output_file, "\n")
fwrite(output, output_file)

cat("Successfully combined", length(file_list), "protein predictions for", nrow(output), "individuals\n")
cat("Output dimensions:", nrow(output), "rows x", ncol(output), "columns\n")