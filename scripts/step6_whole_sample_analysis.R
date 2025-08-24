# ============================================================================
# Step 6: Whole Sample pQTL Analysis
# ============================================================================
# Performs protein quantitative trait loci (pQTL) analysis using the entire 
# sample (no cross-validation). Applies the best-performing statistical method 
# identified in Step 5b to train final models on all available data for 
# subsequent protein prediction and PWAS analysis.
#
# Prerequisites: Must run AFTER Step 5b

# Load packages
suppressMessages({
    library(tidyverse)
    library(plink2R)
    library(pecotmr)
    library(data.table)
    library(future)
    source("./utilities/pqtl_weights.R")
    source("./utilities/pqtl_functions.R")
    source("./utilities/timing_function.R")
})

# ============================================================================
# CONFIGURATION - MODIFY THESE PATHS FOR YOUR ENVIRONMENT
# ============================================================================

# ANALYSIS PARAMETERS
protein_name <- "a1bg"         # Protein name (e.g., "a1bg")
chr_num <- "19"                # Chromosome number (e.g., "19", "01")
min_q_threshold <- 0.2         # FDR threshold for variant inclusion

# INPUT 1: Protein residuals from Step 2
# - Files: {protein_name}_residual_npx.csv
# - Covariate-adjusted NPX residuals for genetic analysis
protein_residuals_dir <- "path_to_protein_residuals"

# INPUT 2: Best methods from Step 5b
# - File: best_method_by_cv_for_protein.csv
# - Contains optimal method for each protein based on CV performance
# - Format: [protein, best_method, best_avg_r2, avg_cv_corr, best_r2_sd]
best_methods_file <- "path_to_best_methods_file"

# INPUT 3: FDR summaries from Step 3
# - LD block and variant selection based on FDR threshold
fdr_summary_dir <- "path_to_step3_fdr_summaries"

# INPUT 4: Genotype data
# - PLINK format files (.bed/.bim/.fam) by chromosome
# - Same genotype data used in previous steps
# - Chromosome-specific analysis for computational efficiency
plink_path <- "path_to_training_genotype_data"

# OUTPUT 1: Final model weights
# - Files: {protein}_{chr}_final_weights.csv
# - Contains trained weights from best method on whole sample
# - Format: [variant_id, weight_value]
# - Used for protein prediction in Steps 7-8
# OUTPUT 2: Model performance
# - Files: {protein}_{chr}_whole_sample_performance.csv
# - Performance metrics on whole sample training
# - Format: [method, r2, correlation, rmse, mae]
output_dir <- "path_to_step6_final_weights"

# END CONFIGURATION
# ============================================================================

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}

# Record start time
time_1 <- Sys.time()

# Load protein residuals
cat("Loading protein residuals for", protein_name, "...\n")
protein_residual_file <- file.path(protein_residuals_dir, paste0(protein_name, "_residual_npx.csv"))

if (!file.exists(protein_residual_file)) {
    stop("Protein residual file not found: ", protein_residual_file)
}

protein_levels <- fread(protein_residual_file)

# Load best methods from Step 5b
cat("Loading best method selection for", protein_name, "...\n")
if (!file.exists(best_methods_file)) {
    stop("Best methods file not found: ", best_methods_file)
}

selected_method_df <- fread(best_methods_file)

selected_method_protein <- selected_method_df %>%
    filter(protein == protein_name)

# Check the number of rows returned
if (nrow(selected_method_protein) == 0) {
    stop(paste("Error: No match found for protein:", protein_name))
} else if (nrow(selected_method_protein) > 1) {
    stop(paste("Error: Multiple rows found for protein:", protein_name))
} else {
    selected_method <- selected_method_protein |>
        pull(best_method)
}

# Read in individual data
protein_levels <- protein_levels |>
    rename_with(~"participant_id", .cols = 1) |>
    rename_with(~"npx", .cols = 2) |>
    arrange(participant_id) |>
    drop_na()

# Read in region summary for FDR flatness (Step 3 output)
fdr_summary_file <- file.path(fdr_summary_dir, paste0(protein_name, "_region_stats_summary.csv"))
if (!file.exists(fdr_summary_file)) {
    stop("FDR summary file not found: ", fdr_summary_file)
}
fdr_summary <- fread(fdr_summary_file)

# List files in region using configurable plink_path
chr_plink_dir <- file.path(plink_path, chr_num)
if (!dir.exists(chr_plink_dir)) {
    stop("Chromosome directory not found: ", chr_plink_dir)
}

region_list <- list.files(
    path = chr_plink_dir,
    pattern = "\\.bim$"
)

if (length(region_list) == 0) {
    stop("No .bim files found in directory: ", chr_plink_dir)
}

## Start loop
for (region in region_list) {
    ## Check for FDR-adjusted p-value < 0.2
    # Determine region start/end

    region_id <- str_extract(region, "^[:graph:]+(?=_un)")

    region_signal <- fdr_summary[region_name == region_id, ][[paste0("signal_FDR", min_q_threshold)]]

    # Formatting output
    print(paste0("========================= Region ", region_id, "=========================="))
    print(paste0("Timestamp: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
    min_p_for_region <- fdr_summary[region_name == region_id, min_pval]
    print(paste0("Smallest p-value: ", min_p_for_region))
    # Note: a couple regions close to the center of the chromosome might have no summary stats
    #       and thus return NA for min_p and min_q. We designate these empty regions 'flat'
    output_p_q_file <- paste0(output_dir, "/", protein_name, "_", region_id, "_minp_values.csv")
    output_p_q <- data.frame(protein_name, region_id, min_p_for_region, region_signal)
    write.table(output_p_q, output_p_q_file, sep = ",", row.names = FALSE)

    if (region_signal == 0 | is.na(region_signal)) {
        print("Flat: TRUE")
    } else if (region_signal == 1) {
        # Print p-value
        print("Flat: FALSE")

        # 1. Import the genotype data
        import_start <- Sys.time()
        region_path <- file.path(chr_plink_dir, tools::file_path_sans_ext(region))
        plink2r_data <- read_plink(region_path)
        gen_matrix <- plink2r_data$bed
        rownames(gen_matrix) <- str_extract(rownames(gen_matrix), "(?<=:)(\\d+)")
        rm(plink2r_data)
        import_end <- Sys.time()
        import_duration <- import_end - import_start
        cat("Step 1. Import plink2 data takes:", format(import_duration, digits = 3), "\n")

        # 2. Match phenotype, genotype, and covariates by IID
        merge_start <- Sys.time()
        common_ids <- intersect(protein_levels$participant_id, rownames(gen_matrix))

        gen_matrix <- gen_matrix[common_ids, , drop = FALSE]
        protein_levels_matched <- protein_levels[protein_levels$participant_id %in% common_ids, ]
        merge_end <- Sys.time()
        merge_duration <- merge_end - merge_start
        # Print durations and start analysis info
        cat("Step 2. Merge pheno/geno/covariates data takes:", format(merge_duration, digits = 3), "\n")
        pqtl_analysis_start <- Sys.time()

        # 3.a Prepare data for analysis
        df_analysis <- cbind(protein_levels, gen_matrix) |>
            janitor::clean_names() # this handles duplicated rsids from multiallelic variants

        # 3.b Assign the selected method to `weight_methods`

        if (selected_method == "bayes_r_weights") {
            weight_methods <- list(bayes_r_weights = list(nit = 1000, nburn = 200, nthin = 5))
        } else if (selected_method == "lasso_weights") {
            weight_methods <- list(lasso_weights = list())
        } else if (selected_method == "enet_weights") {
            weight_methods <- list(enet_weights = list())
        } else if (selected_method == "susie_weights"){
            weight_methods <- list(susie_weights = list())
        } else {
            stop(paste("Error: Invalid method selected:", selected_method))
        }

        ## Update weights methods if only one variant
        if (ncol(gen_matrix) == 1) {
            weight_methods <- list(lm_weights = list())
        }

        post_weights <- list()

        # 3.c Run the selected method on all data (rather than train data)
        post_weights[[1]] <- df_analysis |>
            select(-(1:2)) |>
            as.matrix() |>
            pqtl_weights(
                y = df_analysis$npx,
                X = _,
                weight_methods = weight_methods
            )
        # Convert the matrix to a data frame with appropriate column names
            # Create the data frame based on the new logic for post_weights
            if (!is.null(rownames(post_weights[[1]]))) {
                post_weights_df <- data.frame(
                    variant_id = rownames(post_weights[[1]]),      # Get the row names as variant IDs
                    weights = post_weights[[1]],                    # Use the matrix as weights
                    stringsAsFactors = FALSE                         # Avoid converting strings to factors
                )
            } else if (length(post_weights[[1]]) == length(colnames(df_analysis)[-c(1, 2)])) {
                post_weights_df <- data.frame(
                    variant_id = colnames(df_analysis)[-c(1, 2)],  # Get the column names as variant IDs
                    weights = post_weights[[1]],                    # Use the matrix as weights
                    stringsAsFactors = FALSE                         # Avoid converting strings to factors
                )
            } else if(!all(apply(gen_matrix, 2, \(col) sd(col, na.rm = T) != 0))) {
                post_weights_df <- data.frame(
                    variant_id = colnames(gen_matrix)[apply(gen_matrix, 2, \(col) sd(col, na.rm = T) != 0)],
                    weights = post_weights[[1]],
                    stringsAsFactors = FALSE
                ) # handles the case where variants were dropped due to 0 sd
            } else {
                # Set everything to NA if no conditions are met
                post_weights_df <- data.frame(
                    variant_id = NA,      # Assign NA for variant IDs
                    weights = NA,         # Assign NA for weights
                    stringsAsFactors = FALSE                         # Avoid converting strings to factors
                )
            }
        # Write to a CSV file with specified column names
        write.table(
          post_weights_df,
          file = paste0(
            output_dir, "/", protein_name,
            "_", region_id,
            "_posterior_weights.csv"
          ),
          sep = ",",
          row.names = FALSE,                          # Do not include row numbers
          col.names = TRUE,                           # Include column names
          quote = TRUE                                # Quote character strings
        )

        pqtl_analysis_end <- Sys.time()
        pqtl_analysis_duration <- pqtl_analysis_end - pqtl_analysis_start
        cat("3. pqtl anlaysis for ", length(common_ids), " individuals takes", format(pqtl_analysis_duration, digits = 3), "\n")

        # 4. Run SuSiE using susie_wrapper
        susie_analysis_start <- Sys.time()
        X <- df_analysis %>%
          select(-participant_id, -npx) |>
          as.matrix() |>
          clean_matrix_na() # mode imputation
        susie_res <- pecotmr::susie_wrapper(X = X, y = df_analysis$npx)
        saveRDS(list(susie_res = susie_res), file = paste0(output_dir, "/", protein_name, "_", region_id, "_susie_res.rds"))
        susie_analysis_end <- Sys.time()
        susie_analysis_duration <- susie_analysis_end - susie_analysis_start
        cat("4. susie anlaysis for ", length(df_analysis$npx), " individuals takes", format(susie_analysis_duration, digits = 3), "\n")

        # 5. Free up memory
        rm(gen_matrix, post_weights, protein_levels_matched, df_analysis, X, susie_res)
    } else {
        print(paste0("WARNING: Region ", region_id, " not found in FDR summary."))
    }
}

# Output duration of chromosome analysis
time_2 <- Sys.time()
print(paste0("Running this script takes ", hms_span(time_1, time_2), "."))
