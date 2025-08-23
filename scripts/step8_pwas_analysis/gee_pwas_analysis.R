# ============================================================================
# Step 8: Proteome-Wide Association Study (PWAS) Analysis
# ============================================================================
# 
# Purpose:
# This script performs proteome-wide association studies (PWAS) to test associations
# between predicted protein levels from Step 7 and complex traits/diseases. It uses
# generalized estimating equations (GEE) to account for population structure and
# relatedness in the UK Biobank cohort.
#
# Input:
# 1. Combined proteome predictions from Step 7b:
#    - File: combined_predicted_proteome_all_white_british.csv
#    - Contains predicted NPX values for all proteins across all individuals
#    - Format: [IID, protein1, protein2, ..., proteinN]
#
# 2. Phenotype data:
#    - CSV file with participant phenotypes/traits of interest
#    - Required columns: participant ID, phenotype values, covariates
#    - Format: [IID, phenotype, age, sex, bmi, other_covariates]
#
# 3. Population structure files:
#    - Kinship matrix: for accounting for relatedness
#    - Principal components: for population stratification control
#    - Format: standard PLINK/genetic analysis formats
#
# Output:
# 1. PWAS results: {trait}_pwas_results.csv
#    - Association test results for all proteins with the specified trait
#    - Format: [protein, estimate, std.error, p.value, rho]
#    - estimate: effect size (beta coefficient)
#    - std.error: standard error of the estimate
#    - p.value: statistical significance
#    - rho: correlation parameter from GEE model
#
# Usage Examples:
# For continuous trait: Set pheno_col <- "adj_hdl", analysis_family <- "gaussian"
# For binary trait: Set pheno_col <- "disease_status", analysis_family <- "binomial"
#
# Prerequisites:
# - Must run AFTER Step 7b (combine_npx_files.R)
# - Combined proteome prediction file must exist
# - Phenotype and population structure files must be properly formatted
# ============================================================================

# Load required packages
suppressMessages({
    library(tidyverse)
    library(data.table)
    library(broom.mixed)
    library(igraph)
    library(geepack)
})

# ============================================================================
# CONFIGURATION - MODIFY THESE PATHS FOR YOUR ENVIRONMENT
# ============================================================================

# Phenotype and analysis parameters
pheno_file <- "path_to_phenotype_data.csv"     # Phenotype data file
pheno_col <- "adj_hdl"                         # Phenotype column name (e.g., "adj_hdl", "disease_status")
analysis_family <- "gaussian"                  # Statistical family: "gaussian" (continuous) or "binomial" (binary)

# Input files - modify for your environment
proteome_predictions_file <- "path_to_combined_proteome_predictions.csv"  # Step 7b output
kin_file <- "path_to_kinship_matrix.kin0"      # Kinship matrix file for relatedness correction
pc_file <- "path_to_principal_components.txt"  # Principal components file for population stratification

# Output file
output_file <- "path_to_pwas_results.csv"      # Output file for PWAS results

# ============================================================================
# SCRIPT EXECUTION
# ============================================================================

cat("Starting PWAS analysis...\n")
cat("Phenotype file:", pheno_file, "\n")
cat("Phenotype column:", pheno_col, "\n")
cat("Analysis family:", analysis_family, "\n")
cat("Output file:", output_file, "\n")

# Define Functions
protein_analysis_gee <- function(protein_name){    
    model <- df_pheno |>
        geeglm(as.formula(paste0(pheno_col, " ~ age*sex + bmi + ", protein_name, " + ", paste0("pc", 1:20, collapse = " + "))),
                data = _, id = fam, corstr = "exchangeable", family = analysis_family)
    return(model)
}

clean_analysis <- function(model){
    model |>
        tidy() |>
        filter(term %in% protein_list) |>
        select(estimate, std.error, p.value) |>
        mutate(rho = model$geese$alpha)
}

# Combine functions to improve memory - no storing model
protein_results <- function(protein_name){
    model <- protein_analysis_gee(protein_name)
    results <- clean_analysis(model)
    return(results)
}

# Load Data
cat("Loading proteome predictions...\n")
if (!file.exists(proteome_predictions_file)) {
    stop("Proteome predictions file not found: ", proteome_predictions_file)
}

df_protein <- fread(proteome_predictions_file) |>
    janitor::clean_names()

protein_list <- colnames(df_protein)[-1]

df_analysis <- tibble(protein = protein_list) 

## Kinship data
df_kin <- fread(kin_file) |>
    janitor::clean_names()

## Stroke data
df_pheno <- fread(pheno_file) |>
    janitor::clean_names() |>
    drop_na(all_of(pheno_col))

## PC data
df_pc <- fread(pc_file) |>
    janitor::clean_names()

# Combine stroke data
## Identify kinship
relation_graph <- df_kin |>
    select(iid1, iid2) |>
    graph_from_data_frame(directed = FALSE)

relation_components <- components(relation_graph)

df_fam <- tibble(
    iid = names(relation_components$membership),
    fam = relation_components$membership) |>
    mutate(iid = as.integer(iid))

## Add to stroke data
df_pheno <- left_join(
    df_pheno,
    df_fam, 
    by = join_by(iid)
    ) |>
    mutate(fam = ifelse(is.na(fam), iid, fam)) |>
    left_join(df_pc, by = join_by(iid))

## Add protein data to stroke data
### Inner join effectively filters the stroke dataset to only White British
df_pheno <- inner_join(df_pheno, df_protein, by = join_by(iid))
rm(df_protein, df_fam, df_kin)

print(paste0( "Data imported and processed. Beginning PWAS for ", pheno_col, " with ", nrow(df_pheno), " individuals."))

# Run PWAS with GEE
df_analysis <- df_analysis |>
    mutate(results = map(protein, protein_results)) |>
    unnest(results)

# Output results
df_analysis |>
    janitor::clean_names() |>
    mutate(q_value = p.adjust(p_value, method = "BH")) |>
    arrange(q_value) |>
    saveRDS(paste0(output_file, ".rds"))