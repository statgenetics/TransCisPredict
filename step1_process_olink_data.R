# ============================================================================
# Step 1: Process OLINK Proteomics Data
# ============================================================================
# Processes raw OLINK proteomics data from UK Biobank Pharma Proteomics Project.
# Takes raw batch files containing all proteins and splits them into individual 
# protein files for downstream analysis.

# Load packages
suppressMessages({
    library(tidyverse)
    library(data.table)
    source("./utilities/timing_function.R")
})

# ============================================================================
# CONFIGURATION - MODIFY THESE PATHS FOR YOUR ENVIRONMENT
# ============================================================================

# INPUT: Raw OLINK batch files containing all proteins from UKB-PPP
# - Each batch file contains: participant IDs (eid) + protein expression values (NPX)
# - Multiple batch files numbered sequentially (e.g., 0001.csv, 0002.csv, etc.)
# - Input file format: CSV with columns [eid, protein1, protein2, ..., proteinN]
input_path <- "path_where_you_store_NPX_data_all_proteins"

# Input file naming pattern - use {num} as placeholder for sequential number
# Will be formatted as 4-digit: 0001, 0002, etc.
input_file_template <- "appXXXXX_your_application_ID_olink_instance_0_{num}.csv"

# OUTPUT: Individual protein files, one per protein
# - Each file contains: [eid, protein_expression_value] 
# - Files named: {protein_name}_npx_instance_0_ukb_ppp.csv
# - Missing values removed, participant IDs sorted
output_dir <- "path_where_you_save_each_protein_file"

# END CONFIGURATION
# ============================================================================

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
cat("File template:", input_file_template, "\n")

# Automatically detect available batch files by listing files in input directory
# Extract the pattern from the template to match existing files
file_pattern <- gsub("\\{num\\}", "[0-9]{4}", input_file_template)
available_files <- list.files(input_path, pattern = file_pattern, full.names = FALSE)

if (length(available_files) == 0) {
    stop("No files found matching pattern '", file_pattern, "' in directory: ", input_path)
}

cat("Found", length(available_files), "batch files to process\n")

# Loop through all available batch files
for (file_name in available_files) {
    
    # Extract the number from the filename for logging
    num_match <- regmatches(file_name, regexpr("[0-9]{4}", file_name))
    
    # Construct full input file path
    input_file <- file.path(input_path, file_name)
    
    # Check if file exists before processing
    if (!file.exists(input_file)) {
        cat("Warning: File not found:", input_file, "- skipping\n")
        next
    }
    
    cat("Processing batch file:", num_match, "(", file_name, ")\n")
    
    # Load batch data
    df_protein <- fread(input_file)
    
    # Get list of protein columns (all columns except participant ID)
    protein_list <- df_protein %>%
        select(-eid) %>%
        names()
    
    cat("Found", length(protein_list), "proteins in batch", num_match, "\n")
    
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
    
    cat("Completed batch", num_match, "\n")
}

# Calculate and display processing time
time_2 <- Sys.time()
processing_time <- hms_span(time_1, time_2)

cat("\n============================================================================\n")
cat("OLINK data processing completed successfully!\n")
cat("Total processing time:", processing_time, "\n")
cat("Output files saved to:", output_dir, "\n")
cat("============================================================================\n")

