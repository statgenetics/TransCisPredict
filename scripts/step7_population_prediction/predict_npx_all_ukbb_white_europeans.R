# ============================================================================
# Step 7a: Population-Level Protein Prediction
# ============================================================================
# 
# Purpose:
# This script generates protein expression predictions for the entire UK Biobank
# white European cohort using final model weights from Step 6. It applies trained
# pQTL models to genotype data to predict protein levels for all individuals,
# creating the foundational dataset for proteome-wide association studies (PWAS).
#
# Input:
# 1. Final model weights from Step 6:
#    - Files: {protein}_{chr}_posterior_weights.csv 
#    - Contains trained weights from best method on whole sample
#    - Located in Step 6 output directory
#
# 2. Population genotype data:
#    - PLINK format files (.bed/.bim/.fam) for UK Biobank white Europeans
#    - Complete genotype dataset for target prediction population
#    - Organized by LD blocks matching training data structure
#
# Output:
# 1. Predicted protein levels: {protein}_predicted_npx_all_white_europeans.csv
#    - Contains predicted NPX values for all individuals in target population
#    - Format: [participant_id, predicted_npx]
#    - Used as input for PWAS analysis in Step 8
#
# Usage Example (for a1bg protein):
# Rscript predict_npx_all_ukbb_white_europeans.R a1bg path_to_step6_weights path_to_output
#
# Command line arguments:
# Argument 1: protein_name - Name of protein to predict (e.g., "a1bg", "APOE")
# Argument 2: input_dir - Directory containing final weights from Step 6
# Argument 3: output_dir - Directory to write population predictions
#
# Prerequisites:
# - Must run AFTER Step 6 (whole sample pQTL analysis)
# - Final model weights must exist for the specified protein
# - Population genotype data must be available and formatted consistently
#
# Note: This step requires substantial memory (50-100GB) for large-scale prediction
# ============================================================================

# Clear environment and record start time
rm(list = ls())
time_1 <- Sys.time()

## CONFIGURATION - MODIFY THESE PATHS FOR YOUR ENVIRONMENT

# Command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
    stop("Usage: Rscript predict_npx_all_ukbb_white_europeans.R <protein_name> <weights_dir> <output_dir>")
}

protein_name <- args[1]    # Protein name (e.g., "a1bg", "APOE")
input_dir <- args[2]       # Directory with final weights from Step 6
output_dir <- args[3]      # Output directory for predictions

# Start time to later measure total run time
time_1 <- Sys.time()

# Input paths - modify for your environment
population_genotype_path <- "path_to_population_genotype_data"  # UK Biobank white Europeans genotypes

# Load packages
suppressMessages({
    library(tidyverse)
    library(data.table)
    library(plink2R)
    library(rsample)
    source("../utilities/pqtl_functions.R")
    source("../utilities/timing_function.R")
})

# Create output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("Starting protein prediction for", protein_name, "...\n")
cat("Weights directory:", input_dir, "\n")
cat("Output directory:", output_dir, "\n")

# Find final weight files for this protein
post_weights_list <- list.files(input_dir, pattern = paste0(protein_name, ".*posterior_weights\\.csv$"))

if (length(post_weights_list) == 0) {
    stop("No weight files found for protein ", protein_name, " in directory: ", input_dir)
}

cat("Found", length(post_weights_list), "weight files for", protein_name, "\n")

# Initialize output - adjust population size as needed
population_size <- 407917  # UK Biobank white Europeans sample size
results <- vector(length = population_size)

## PART 3: REGION CONTRIBUTION
# Loop over the directory of posterior weights
for (file in post_weights_list){

    # Start region
    region <-  str_extract(file, paste0("(?<=", protein_name, "_)[:graph:]+(?=_post)"))
    print(paste0("Starting region ", region))
    
    # Load in and prepare weights
    weights <- fread(paste0(input_dir, file)) |>
        drop_na()
    
    weights_matrix <- as.matrix(weights[, -1])
    rownames(weights_matrix) <- weights$variant_id
    
    # Load in and prepare genotype
    plink2r_data <- read_plink(paste0( 
            "/mnt/vast/hpc/csg/NotBackedUp/rd2972/20240323_UKBB_proteomics/20240909_UKBB_white_Europeans_407917_individuals/genotyped_variants/",
            "UKBB_white_Europeans_407917_individuals_", region))
    
    # Extract genotype matrix and remove plink data
    gen_matrix <- plink2r_data$bed
    rownames(gen_matrix) <- str_extract(rownames(gen_matrix), "(?<=:)(\\d+)")
    iid_vec <- plink2r_data$fam$V1
    fid_vec <- plink2r_data$fam$V2
    rm(plink2r_data)
    
    # QC the columns
    rsid_overlap = intersect(rownames(weights_matrix), colnames(gen_matrix))
    if (length(rsid_overlap) == 0) stop("Misalignment between weights and genotype")
    
    weights_matrix <- weights_matrix[rsid_overlap,,drop=FALSE]
    gen_matrix <- gen_matrix[,rsid_overlap, drop=FALSE]
    
    gen_matrix_imputed <- apply(gen_matrix, 2, function(col) {
    if (any(is.na(col))) {
        col[is.na(col)] <- stat_mode(col[!is.na(col)])
    }
        return(col)
    })
    
    # Add region's contribution to output
    results <- results + (gen_matrix_imputed %*% weights_matrix)
    
    # Remove excess objects
    rm(gen_matrix, gen_matrix_imputed)
}

# PART 4: PRINT OUTPUT
output <- tibble(
    IID = iid_vec,
    FID = fid_vec,
    npx = results[,1])

write_csv(output, file = paste0(output_dir, protein_name, "_npx_ukb_407917_white_europeans.csv"))

# Print time of calculation
time_2 <- Sys.time()
print(paste0("Running this script takes ", hms_span(time_1, time_2), "."))
