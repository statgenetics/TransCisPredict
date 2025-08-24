# ============================================================================
# Step 2: Covariate Regression for Protein Expression
# ============================================================================
# Removes covariate effects from protein expression levels by performing linear 
# regression and extracting residuals. Generates covariate-adjusted protein 
# phenotypes for genetic analysis while preserving genetic signal.
#
# Statistical Model: NPX ~ β₀ + β1×age + β2×sex + β3×age×sex + β4×BMI + Σ(βi×PCi) + ε
# Residuals = NPX - (predicted values from covariates)

# Load packages
suppressMessages({
    library(tidyverse)
    library(data.table)
    library(broom)
    library(janitor)
})

# ============================================================================
# CONFIGURATION - MODIFY THESE PATHS FOR YOUR ENVIRONMENT
# ============================================================================

# INPUT 1: Individual protein files from Step 1
# - Directory containing files: {protein_name}_npx_instance_0_ukb_ppp.csv
# - Each file format: [eid, protein_expression_value]
# - Generated from Step 1 data processing
protein_files_dir <- "path_where_you_save_each_protein_file"

# INPUT 2: Covariate file
# - CSV file containing participant covariates
# - Required columns: IID (participant ID), age, sex, bmi, pc1-pc20
# - Format: [IID, age, sex, bmi, pc1, pc2, ..., pc20, other_covariates]
# - Model will include: age + sex + age×sex + BMI + PC1-20
covariate_file <- "path_to_your_covariate_file.csv"

# OUTPUT 1: Residual protein expression files (key input for Steps 4-6)
# - Files: {protein_name}_residual_npx.csv
# - Format: [eid, npx_residuals]
output_residuals_dir <- "path_to_save_protein_residuals"

# OUTPUT 2: Covariate coefficient files
# - Files: {protein_name}_covariate_coefficients.csv  
# - Format: [term, estimate, std.error, p.value]
output_coefficients_dir <- "path_to_save_covariate_coefficients"

# OUTPUT 3: Model summary file
# - File: protein_covariate_r2_summary.csv
# - Format: [protein, r_squared, n_individuals]
output_summary_file <- "protein_covariate_r2_summary.csv"

# END CONFIGURATION
# ============================================================================

# Create output directories if they don't exist
if (!dir.exists(output_residuals_dir)) {
    dir.create(output_residuals_dir, recursive = TRUE)
}
if (!dir.exists(output_coefficients_dir)) {
    dir.create(output_coefficients_dir, recursive = TRUE)
}

# ============================================================================
# DATA LOADING AND PREPARATION
# ============================================================================

cat("Starting covariate regression analysis...\n")
cat("Protein files directory:", protein_files_dir, "\n")
cat("Covariate file:", covariate_file, "\n")
cat("Output directories:\n")
cat("  - Residuals:", output_residuals_dir, "\n")  
cat("  - Coefficients:", output_coefficients_dir, "\n")

# Get list of protein files from Step 1 output
protein_files <- list.files(protein_files_dir, pattern = "_npx_instance_0_ukb_ppp.csv$", full.names = FALSE)

if (length(protein_files) == 0) {
    stop("No protein files found in directory: ", protein_files_dir)
}

# Extract protein names from filenames
protein_list <- protein_files %>%
    str_extract("^[^_]+(?=_npx_instance_0_ukb_ppp\\.csv$)")

cat("Found", length(protein_list), "protein files to process\n")

# Load covariate data
if (!file.exists(covariate_file)) {
    stop("Covariate file not found: ", covariate_file)
}

covariates_ref <- fread(covariate_file)
cat("Loaded covariate data for", nrow(covariates_ref), "individuals\n")

# Verify required columns exist in covariate file
required_cols <- c("IID", "age", "sex", "bmi")
missing_cols <- setdiff(required_cols, names(covariates_ref))
if (length(missing_cols) > 0) {
    stop("Missing required columns in covariate file: ", paste(missing_cols, collapse = ", "))
}

# Check for principal components (PC1-PC20)
pc_cols <- paste0("pc", 1:20)
available_pcs <- intersect(pc_cols, names(covariates_ref))
cat("Found", length(available_pcs), "principal components in covariate file\n")

if (length(available_pcs) < 10) {
    warning("Fewer than 10 principal components found. Consider including more PCs for population structure control.")
}

# ============================================================================
# COVARIATE REGRESSION PROCESSING
# ============================================================================

# Initialize results tracking
r2_list <- c()
processed_proteins <- c()
failed_proteins <- c()

cat("\nProcessing proteins...\n")

# Loop through each protein
for (i in seq_along(protein_list)) {
    protein_name <- protein_list[i]
    
    if (i %% 100 == 0) {
        cat("Processing protein", i, "of", length(protein_list), ":", protein_name, "\n")
    }
    
    tryCatch({
        # Load protein NPX data
        protein_file_path <- file.path(protein_files_dir, protein_files[i])
        npx <- fread(protein_file_path)
        
        # Check file format
        if (ncol(npx) != 2) {
            stop("Unexpected file format for ", protein_name, ". Expected 2 columns, found ", ncol(npx))
        }
        
        # Standardize column names
        names(npx) <- c("eid", "npx_value")
        
        # Find overlapping individuals between protein and covariate data
        overlap_ids <- intersect(npx$eid, covariates_ref$IID)
        
        if (length(overlap_ids) < 100) {
            warning("Only ", length(overlap_ids), " overlapping individuals for ", protein_name, ". Skipping.")
            failed_proteins <- c(failed_proteins, protein_name)
            next
        }
        
        # Filter to overlapping individuals
        npx_filtered <- npx %>%
            filter(eid %in% overlap_ids) %>%
            arrange(eid)
            
        covariates_filtered <- covariates_ref %>%
            filter(IID %in% overlap_ids) %>%
            arrange(IID)
        
        # Verify alignment between protein and covariate data
        if (!identical(npx_filtered$eid, covariates_filtered$IID)) {
            stop("Misalignment between protein and covariate participant IDs for ", protein_name)
        }
        
        # Prepare regression data frame
        regression_df <- covariates_filtered %>%
            mutate(npx = npx_filtered$npx_value) %>%
            select(-IID) %>%
            clean_names()  # Standardize column names
        
        # Build regression formula with available covariates
        covariate_terms <- c("age", "sex", "bmi")
        
        # Add available principal components  
        available_pcs_clean <- intersect(paste0("pc", 1:20), names(regression_df))
        covariate_terms <- c(covariate_terms, available_pcs_clean)
        
        # Add age*sex interaction if both are available
        if (all(c("age", "sex") %in% names(regression_df))) {
            covariate_terms <- c(covariate_terms, "I(age * sex)")
        }
        
        # Create regression formula
        formula_str <- paste("npx ~", paste(covariate_terms, collapse = " + "))
        
        # Fit linear regression model
        model <- lm(as.formula(formula_str), data = regression_df)
        
        # Extract residuals
        residuals_df <- tibble(
            eid = npx_filtered$eid,
            npx_residuals = residuals(model)
        )
        
        # Save residual protein expression (key output for downstream analysis)
        residuals_output_path <- file.path(output_residuals_dir, paste0(protein_name, "_residual_npx.csv"))
        write_csv(residuals_df, residuals_output_path)
        
        # Save model coefficients
        coefficients_df <- model %>%
            tidy() %>%
            select(term, estimate, std.error, p.value)
            
        coefficients_output_path <- file.path(output_coefficients_dir, paste0(protein_name, "_covariate_coefficients.csv"))
        write_csv(coefficients_df, coefficients_output_path)
        
        # Store R-squared value
        model_summary <- summary(model)
        r2_value <- model_summary$r.squared
        r2_list <- c(r2_list, r2_value)
        processed_proteins <- c(processed_proteins, protein_name)
        
    }, error = function(e) {
        cat("Error processing", protein_name, ":", e$message, "\n")
        failed_proteins <- c(failed_proteins, protein_name)
    })
}

# ============================================================================
# SUMMARY AND OUTPUT
# ============================================================================

cat("\n============================================================================\n")
cat("Covariate regression analysis completed!\n")
cat("Successfully processed:", length(processed_proteins), "proteins\n")
cat("Failed proteins:", length(failed_proteins), "\n")

if (length(failed_proteins) > 0) {
    cat("Failed protein list:\n")
    cat(paste(failed_proteins, collapse = ", "), "\n")
}

# Save R-squared summary
if (length(processed_proteins) > 0) {
    r2_summary <- tibble(
        protein = processed_proteins,
        r_squared = r2_list,
        n_individuals = rep(length(overlap_ids), length(processed_proteins))  # Approximate
    )
    
    write_csv(r2_summary, output_summary_file)
    
    cat("\nR-squared summary statistics:\n")
    cat("Mean R²:", round(mean(r2_list), 4), "\n")
    cat("Median R²:", round(median(r2_list), 4), "\n")
    cat("Range R²:", round(min(r2_list), 4), "to", round(max(r2_list), 4), "\n")
    cat("Summary saved to:", output_summary_file, "\n")
}

cat("\nOutput files:\n")
cat("- Protein residuals (", length(processed_proteins), "files ):", output_residuals_dir, "\n")
cat("- Covariate coefficients (", length(processed_proteins), "files ):", output_coefficients_dir, "\n")
cat("============================================================================\n")

