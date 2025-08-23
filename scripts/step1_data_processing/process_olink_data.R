# ============================================================================
# Step 1: Process OLINK Proteomics Data
# ============================================================================
# 
# Purpose:
# This script processes raw OLINK proteomics data from UK Biobank Pharma 
# Proteomics Project (UKB-PPP). It takes the raw batch files containing all 
# proteins and splits them into individual protein files for downstream analysis.
#
# Input:
# - Raw OLINK batch files containing all proteins from UKB-PPP
# - Each batch file contains: participant IDs (eid) + protein expression values (NPX)
# - Multiple batch files numbered sequentially (e.g., 0001.csv, 0002.csv, etc.)
# - Input file format: CSV with columns [eid, protein1, protein2, ..., proteinN]
#
# Output:
# - Individual protein files, one per protein
# - Each file contains: [eid, protein_expression_value] 
# - Files named: {protein_name}_npx_instance_0_ukb_ppp.csv
# - Missing values removed, participant IDs sorted
#
# Processing:
# - Loops through all batch files
# - Extracts each protein column individually
# - Removes missing values (drop_na)
# - Sorts by participant ID for consistency
# - Saves as individual CSV files for downstream analysis
#
# Usage:
# Set input_path to directory containing numbered batch files
# Set output_dir to directory where individual protein files will be saved
# ============================================================================

# Load packages and functions
suppressMessages({
    library(tidyverse)
    library(data.table)
    source("../utilities/timing_function.R")  # Timing utility function
})

# ============================================================================
# CONFIGURATION - MODIFY THESE PATHS FOR YOUR ENVIRONMENT
# ============================================================================

# Input directory containing raw OLINK batch files
# Files should be named with sequential numbers: 0001.csv, 0002.csv, etc.
input_path <- "path_where_you_store_NPX_data_all_proteins"

# Output directory for individual protein files
output_dir <- "path_where_you_save_each_protein_file"

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}

# ============================================================================
# DATA PROCESSING
# ============================================================================

# Mark start time for processing duration tracking
time_1 <- Sys.time()

cat("Starting OLINK proteomics data processing...\n")
cat("Input directory:", input_path, "\n")
cat("Output directory:", output_dir, "\n")

# Loop through numbered batch files (assuming files are numbered 0001-9999)
# Modify the range (1:12) based on your actual number of batch files
for (num in sprintf("%04d", 1:12)) {
    
    # Construct input file path
    input_file <- file.path(input_path, paste0("app32285_20231207101836_olink_instance_0_", num, ".csv"))
    
    # Check if file exists before processing
    if (!file.exists(input_file)) {
        cat("Warning: File not found:", input_file, "- skipping\n")
        next
    }
    
    cat("Processing batch file:", num, "\n")
    
    # Load batch data
    df_protein <- fread(input_file)
    
    # Get list of protein columns (all columns except participant ID)
    protein_list <- df_protein %>%
        select(-eid) %>%
        names()
    
    cat("Found", length(protein_list), "proteins in batch", num, "\n")
    
    # Process each protein individually
    for (protein in protein_list) {
        
        # Extract protein data with participant IDs
        df_working <- df_protein %>%
            select(eid, all_of(protein)) %>%
            drop_na() %>%                    # Remove missing values
            arrange(eid)                     # Sort by participant ID
        
        # Generate output filename
        output_file <- file.path(output_dir, paste0(protein, "_npx_instance_0_ukb_ppp.csv"))
        
        # Write individual protein file
        write_csv(df_working, output_file)
    }
    
    cat("Completed batch", num, "\n")
}

# Calculate and display processing time
time_2 <- Sys.time()
processing_time <- hms_span(time_1, time_2)

cat("\n============================================================================\n")
cat("OLINK data processing completed successfully!\n")
cat("Total processing time:", processing_time, "\n")
cat("Output files saved to:", output_dir, "\n")
cat("============================================================================\n")

# ============================================================================
# EXPECTED OUTPUT
# ============================================================================
#
# Individual protein files with format:
# - Filename: {protein_name}_npx_instance_0_ukb_ppp.csv
# - Columns: eid, {protein_name}
# - Example content:
#   eid,PROTEIN1
#   1000001,5.2345
#   1000002,4.8901
#   1000003,6.1234
#
# Quality control notes:
# - Missing values removed for each protein
# - Participant IDs sorted consistently
# - Ready for Step 2 covariate regression analysis
# ============================================================================