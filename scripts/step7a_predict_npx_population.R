# ============================================================================
# Step 7a: Population-Level Protein Prediction
# ============================================================================
# Generates protein expression predictions for the entire UK Biobank White British 
# cohort using final model weights from Step 6. Applies trained pQTL models to 
# genotype data to predict protein levels for all individuals, creating the 
# foundational dataset for proteome-wide association studies (PWAS).
#
# Prerequisites: Must run AFTER Step 6
# Note: This step requires substantial memory (50-100GB) for large-scale prediction

# Load packages
suppressMessages({
    library(tidyverse)
    library(data.table)
    library(plink2R)
    library(rsample)
    source("./utilities/pqtl_functions.R")
    source("./utilities/timing_function.R")
})

# ============================================================================
# CONFIGURATION - MODIFY THESE PATHS FOR YOUR ENVIRONMENT
# ============================================================================

# ANALYSIS PARAMETER
protein_name <- "a1bg"     # Protein name (e.g., "a1bg")

# INPUT 1: Final model weights from Step 6
# - Files: {protein}_{chr}_posterior_weights.csv 
# - Contains trained weights from best method on whole sample
input_dir <- "path_to_step6_final_weights"

# INPUT 2: Population genotype data
# - PLINK format files (.bed/.bim/.fam) for UK Biobank White British
# - Complete genotype dataset for target prediction population
# - Organized by LD blocks matching training data structure
# - Base path to population genotype data (without region suffix)
# - Script will append "_{region}" to construct full PLINK file paths
genotype_path <- "path_to_population_genotype_data"

# OUTPUT: Predicted protein levels
# - File: {protein}_predicted_npx_all_white_british.csv
# - Contains predicted NPX values for all individuals in target population
# - Format: [participant_id, predicted_npx]
# - Used as input for PWAS analysis in Step 8
output_dir <- "path_to_step7_population_predictions"

# END CONFIGURATION
# ============================================================================

# Clear environment and record start time
rm(list = ls()[!ls() %in% c("protein_name", "input_dir", "genotype_path", "output_dir")])
time_1 <- Sys.time()

# Create output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("Starting protein prediction for", protein_name, "...\n")
cat("Weights directory:", input_dir, "\n")
cat("Output directory:", output_dir, "\n")
cat("Genotype path:", genotype_path, "\n")

# Find final weight files for this protein
post_weights_list <- list.files(input_dir, pattern = paste0(protein_name, ".*posterior_weights\\.csv$"))

if (length(post_weights_list) == 0) {
    stop("No weight files found for protein ", protein_name, " in directory: ", input_dir)
}

cat("Found", length(post_weights_list), "weight files for", protein_name, "\n")

# Determine population size from .fam file
# Read the first available genotype file to get population size
first_weight_file <- post_weights_list[1]
first_region <- str_extract(first_weight_file, paste0("(?<=", protein_name, "_)[:graph:]+(?=_post)"))
temp_plink_data <- read_plink(paste0(genotype_path, "_", first_region))
population_size <- nrow(temp_plink_data$fam)
rm(temp_plink_data)  # Clean up temporary data

cat("Population size determined from .fam file:", population_size, "\n")

# Initialize output 
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
    plink2r_data <- read_plink(paste0(genotype_path, "_", region))
    
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

write_csv(output, file = paste0(output_dir, protein_name, "_npx_ukb_407917_white_british.csv"))

# Print time of calculation
time_2 <- Sys.time()
print(paste0("Running this script takes ", hms_span(time_1, time_2), "."))
