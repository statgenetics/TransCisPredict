# ============================================================================
# Step 7b: Combine Population Protein Predictions 
# ============================================================================
# 
# Purpose:
# This script combines individual protein prediction files from Step 7a into a
# single comprehensive dataset containing predicted protein levels for all proteins
# across the entire UK Biobank white European cohort. This creates the final
# proteome-wide prediction matrix for PWAS analysis.
#
# Input:
# 1. Individual protein prediction files from Step 7a:
#    - Files: {protein}_predicted_npx_all_white_europeans.csv
#    - Each contains: [participant_id, predicted_npx] for one protein
#    - Located in Step 7a output directory
#
# Output:
# 1. Combined proteome predictions: combined_predicted_proteome_all_white_europeans.csv
#    - Contains predicted NPX values for ALL proteins and ALL individuals
#    - Format: [IID, protein1, protein2, protein3, ..., proteinN]
#    - One row per individual, one column per protein
#    - Used as input for Step 8 PWAS analysis
#
# Usage Example:
# Rscript combine_npx_files.R path_to_step7a_predictions path_to_combined_output.csv
#
# Command line arguments:
# Argument 1: input_dir - Directory containing individual protein prediction files from Step 7a
# Argument 2: output_file - Full path to combined output file (including filename)
#
# Prerequisites:
# - Must run AFTER Step 7a (predict_npx_all_ukbb_white_europeans.R)
# - All protein prediction files must have consistent participant ID formatting
# - Files should contain the same number of individuals for proper merging
# ============================================================================

# Load required packages
suppressMessages({
    library(tidyverse)
    library(data.table)
})

## CONFIGURATION - MODIFY THESE PATHS FOR YOUR ENVIRONMENT

# Command line arguments
args = commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
    stop("Usage: Rscript combine_npx_files.R <input_directory> <output_file>")
}

input_dir = args[1]     # Directory with individual protein predictions from Step 7a
output_file = args[2]   # Full path to combined output file

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