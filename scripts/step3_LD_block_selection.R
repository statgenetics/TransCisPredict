# ============================================================================
# Step 3: LD Block Selection Based on FDR Thresholds
# ============================================================================
# Performs marginal genetic association analysis across all LD blocks to identify 
# regions with significant genetic signal. Tests each variant within LD blocks and 
# classifies blocks based on FDR-corrected significance thresholds. This filtering 
# step focuses downstream analysis on informative genomic regions.
#
# Statistical Model (per SNP): NPX_residuals ~ age + sex + age×sex + BMI + PC1-20 + SNP + ε

# Load packages
suppressMessages({
    library(data.table)
    library(tidyverse)
    library(plink2R)
    library(broom)
    library(R.utils)
    source("./utilities/pqtl_functions.R")
    source("./utilities/pqtl_weights.R")
    source("./utilities/timing_function.R")
})

# ============================================================================
# CONFIGURATION - MODIFY THESE PATHS FOR YOUR ENVIRONMENT
# ============================================================================

# ANALYSIS PARAMETER: Protein to analyze
# - Name of protein to analyze (e.g., "a1bg")
protein_name <- "protein_name"

# INPUT 1: Protein residuals from Step 2
# - Files: {protein_name}_residual_npx.csv  
# - Format: [eid, npx_residuals]
# - Covariate-adjusted protein phenotypes
protein_residuals_dir <- "path_to_save_protein_residuals"

# INPUT 2: Covariate file (same as Step 2)
# - CSV file with: IID, age, sex, bmi, pc1-pc20
# - Used for marginal association testing
covariate_file <- "path_to_your_covariate_file.csv"

# INPUT 3: Genotype data by LD blocks
# - PLINK format files (.bed/.bim/.fam) split by genomic regions
# - Directory structure: chromosome/LD_block files
# - Genotyped variants only (no imputation)
genotype_base_path <- "path_to_genotype_data_by_LD_blocks"

# OUTPUT: Region summary file
# - File: {protein_name}_region_stats_summary.csv
# - Format: [chr, region_start, region_end, region_name, n_unfiltered_variants, 
#           n_filtered_variants, min_pval, signal_FDR0.1, signal_FDR0.2]
# - Used in Steps 4-5 to focus analysis on informative regions
output_dir <- "path_to_save_LD_block_results"

# ANALYSIS PARAMETERS: FDR thresholds for region classification
fdr_threshold_strict <- 0.1   # Stricter threshold
fdr_threshold_liberal <- 0.2  # More liberal threshold

# END CONFIGURATION
# ============================================================================

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}

cat("============================================================================\n")
cat("LD Block Selection Analysis\n")
cat("============================================================================\n")
cat("Protein:", protein_name, "\n")
cat("Residuals directory:", protein_residuals_dir, "\n")
cat("Genotype base path:", genotype_base_path, "\n")
cat("Output directory:", output_dir, "\n")
cat("FDR thresholds:", fdr_threshold_strict, "and", fdr_threshold_liberal, "\n")

# ============================================================================
# FUNCTIONS
# ============================================================================

# Function to perform marginal association analysis for a single SNP
marginal_analysis <- function(snp_name, phenotype_data){
    
    # Build regression formula including all covariates
    formula_str <- paste0("npx ~ age + sex + age: sex + bmi + ", 
                         paste0("pc", 1:20, collapse = " + "), 
                         " + ", snp_name)
    
    # Fit linear model
    tryCatch({
        model <- lm(as.formula(formula_str), data = phenotype_data)
        
        # Extract p-value for the SNP term
        p_val <- model %>%
            tidy() %>%
            filter(term == snp_name) %>%
            pull(p.value)
        
        return(if(length(p_val) == 0) NA else p_val)
        
    }, error = function(e) {
        return(NA)
    })
}

# ============================================================================
# DATA LOADING
# ============================================================================

cat("\nLoading input data...\n")

# Load protein residuals from Step 2
protein_residual_file <- file.path(protein_residuals_dir, paste0(protein_name, "_residual_npx.csv"))

if (!file.exists(protein_residual_file)) {
    stop("Protein residual file not found: ", protein_residual_file)
}

protein_levels <- fread(protein_residual_file) %>%
    rename(IID = eid, npx = npx_residuals) %>%
    arrange(IID) %>%
    drop_na()

cat("Loaded protein residuals for", nrow(protein_levels), "individuals\n")

# Load covariate data (same as used in Step 2)
if (!file.exists(covariate_file)) {
    stop("Covariate file not found: ", covariate_file)
}

covariates <- fread(covariate_file)
cat("Loaded covariate data for", nrow(covariates), "individuals\n")

# Merge protein and covariate data
df_covariates <- inner_join(protein_levels, covariates, by = join_by(IID))
cat("Combined data:", nrow(df_covariates), "individuals with both protein and covariate data\n")

# Verify required columns
required_cols <- c("IID", "npx", "age", "sex", "bmi")
pc_cols <- paste0("pc", 1:20)
available_pcs <- intersect(pc_cols, names(df_covariates))

missing_cols <- setdiff(c(required_cols, available_pcs), names(df_covariates))
if (length(missing_cols) > 0) {
    warning("Missing columns: ", paste(missing_cols, collapse = ", "))
}

cat("Available principal components:", length(available_pcs), "\n")

# ============================================================================
# ANALYSIS SETUP
# ============================================================================

# Initialize output file with headers
output_file <- file.path(output_dir, paste0(protein_name, "_region_stats_summary.csv"))

# Create header for output file
output_header <- data.frame(
    chr = integer(),
    region_start = integer(),
    region_end = integer(),
    region_name = character(),
    n_unfiltered_variants = integer(),
    n_filtered_variants = integer(),
    min_pval = numeric(),
    signal_FDR0.1 = integer(),
    signal_FDR0.2 = integer()
)

write.csv(output_header, file = output_file, row.names = FALSE, quote = FALSE)
cat("\nInitialized output file:", output_file, "\n")

cat("\nStarting LD block analysis across chromosomes...\n")

total_regions <- 0
regions_with_signal_01 <- 0
regions_with_signal_02 <- 0

# ============================================================================
# CHROMOSOME ITERATION
# ============================================================================

# Process each chromosome
for (chr_num in sprintf("%02d", 1:22)) {
    chr_start_time <- Sys.time()
    
    cat("Processing chromosome", chr_num, "...\n")
    
    # Get LD block directory for this chromosome
    chr_dir <- file.path(genotype_base_path, chr_num)
    
    if (!dir.exists(chr_dir)) {
        cat("Warning: Chromosome directory not found:", chr_dir, "- skipping\n")
        next
    }
    
    # Get list of LD regions (.bim files indicate PLINK filesets)
    region_files <- list.files(chr_dir, pattern = "\\.bim$")
    
    if (length(region_files) == 0) {
        cat("Warning: No .bim files found in", chr_dir, "- skipping\n")
        next
    }
    
    cat("Found", length(region_files), "LD blocks in chromosome", chr_num, "\n")
    
    # Process each LD region
    for (region_file in region_files) {
        
        region_basename <- tools::file_path_sans_ext(region_file)
        region_path <- file.path(chr_dir, region_basename)
        
        tryCatch({
            # Load genotype data using plink2R
            plink_data <- read_plink(region_path)
            gen_matrix <- plink_data$bed
            
            # Extract participant IDs and clean SNP names
            rownames(gen_matrix) <- str_extract(rownames(gen_matrix), "(?<=:)(\\d+)")
            colnames(gen_matrix) <- janitor::make_clean_names(colnames(gen_matrix))
            
            # Get region coordinates
            n_unfiltered_variants <- ncol(gen_matrix)
            region_start <- min(plink_data$bim$V4)
            region_end <- max(plink_data$bim$V4)
            region_name <- str_extract(region_file, "^[^_]+(?=_unrel|\\.)") # Extract clean region name
            
            # Filter to common individuals
            common_ids <- intersect(df_covariates$IID, as.numeric(rownames(gen_matrix)))
            
            if (length(common_ids) < 100) {
                cat("Warning: Only", length(common_ids), "overlapping individuals for", region_name, "- skipping\n")
                next
            }
            
            # Subset data to common individuals
            df_protein <- df_covariates %>% filter(IID %in% common_ids)
            gen_matrix <- gen_matrix[as.character(common_ids), , drop = FALSE]
            
            # Quality control: impute missing values and filter constant variants
            gen_imputed <- apply(gen_matrix, 2, function(col) {
                if (any(is.na(col))) {
                    col[is.na(col)] <- stat_mode(col[!is.na(col)])
                }
                return(col)
            })
            
            # Filter out variants with zero variance
            valid_variants <- apply(gen_imputed, 2, function(col) sd(col, na.rm = TRUE) != 0)
            gen_filtered <- gen_imputed[, valid_variants, drop = FALSE]
            n_filtered_variants <- ncol(gen_filtered)
            
            # Skip if no valid variants remain
            if (n_filtered_variants == 0) {
                cat("Warning: No valid variants in", region_name, "after QC - skipping\n")
                
                # Still record the region with NA values
                region_result <- data.frame(
                    chr = as.numeric(chr_num),
                    region_start = region_start,
                    region_end = region_end,
                    region_name = region_name,
                    n_unfiltered_variants = n_unfiltered_variants,
                    n_filtered_variants = 0,
                    min_pval = NA,
                    signal_FDR0.1 = 0,
                    signal_FDR0.2 = 0
                )
                
                write.table(region_result, file = output_file, row.names = FALSE, 
                           col.names = FALSE, append = TRUE, sep = ",", quote = FALSE)
                next
            }
            
            # Combine genotype and phenotype data
            gen_df <- gen_filtered %>%
                as_tibble() %>%
                mutate(IID = as.numeric(rownames(gen_filtered)))
            
            analysis_df <- df_protein %>%
                inner_join(gen_df, by = "IID")
            
            # Perform marginal association analysis for each SNP
            snp_names <- colnames(gen_filtered)
            
            cat("  Analyzing", length(snp_names), "variants in", region_name, "\n")
            
            # Run association tests
            p_values <- map_dbl(snp_names, ~marginal_analysis(.x, analysis_df))
            p_values <- p_values[!is.na(p_values)]  # Remove failed tests
            
            # Calculate statistics if we have valid p-values
            if (length(p_values) > 0) {
                min_p <- min(p_values, na.rm = TRUE)
                
                # Apply Benjamini-Hochberg FDR correction
                q_values <- p.adjust(p_values, method = "BH")
                min_q <- min(q_values, na.rm = TRUE)
                
                # Classify region based on FDR thresholds
                signal_FDR01 <- ifelse(min_q <= fdr_threshold_strict, 1, 0)
                signal_FDR02 <- ifelse(min_q <= fdr_threshold_liberal, 1, 0)
                
            } else {
                min_p <- NA
                signal_FDR01 <- 0
                signal_FDR02 <- 0
            }
            
            # Record results
            region_result <- data.frame(
                chr = as.numeric(chr_num),
                region_start = region_start,
                region_end = region_end,
                region_name = region_name,
                n_unfiltered_variants = n_unfiltered_variants,
                n_filtered_variants = n_filtered_variants,
                min_pval = min_p,
                signal_FDR0.1 = signal_FDR01,
                signal_FDR0.2 = signal_FDR02
            )
            
            # Append to output file
            write.table(region_result, file = output_file, row.names = FALSE,
                       col.names = FALSE, append = TRUE, sep = ",", quote = FALSE)
            
            # Update counters
            total_regions <- total_regions + 1
            if (signal_FDR01 == 1) regions_with_signal_01 <- regions_with_signal_01 + 1
            if (signal_FDR02 == 1) regions_with_signal_02 <- regions_with_signal_02 + 1
            
        }, error = function(e) {
            cat("Error processing region", region_file, ":", e$message, "\n")
        })
    }
    
    chr_end_time <- Sys.time()
    chr_duration <- hms_span(chr_start_time, chr_end_time)
    cat("Completed chromosome", chr_num, "- Duration:", chr_duration, "\n")
}

# ============================================================================
# SUMMARY AND OUTPUT
# ============================================================================

cat("\n============================================================================\n")
cat("LD Block Selection Analysis Completed!\n")
cat("============================================================================\n")
cat("Protein:", protein_name, "\n")
cat("Total LD blocks analyzed:", total_regions, "\n")
cat("Blocks with signal (FDR < 0.1):", regions_with_signal_01, 
    sprintf(" (%.1f%%)", 100 * regions_with_signal_01 / max(total_regions, 1)), "\n")
cat("Blocks with signal (FDR < 0.2):", regions_with_signal_02, 
    sprintf(" (%.1f%%)", 100 * regions_with_signal_02 / max(total_regions, 1)), "\n")
cat("Output file:", output_file, "\n")
cat("============================================================================\n")

