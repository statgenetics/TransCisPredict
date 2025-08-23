#' QC pQTL Data
#'
#' This function takes a dataframe of summary statistics and an LD matrix then proccesses it 
#' with the pecotmr pipeline. It constructs the necessary data structures and passes them to 
#' the QC functions from pecotmr.
#'
#' @param sumstats A dataframe of summary statistics including either variant_id or chrom, pos, A2, A1 
#'        columns.
#' @param LD An LD matrix with column names in chrom:pos:A2:A1 format.
#' @import tibble dplyr
#' @export
pqtl_qc <- function(sumstats, LD, ...){
    # align row/col names
    rownames(LD) <- colnames(LD)
    
    # convert LD to LD_data
    ld_variant_id = colnames(LD)
    
    ref_panel <- do.call(rbind, lapply(strsplit(ld_variant_id, ":"), function(x) {
        data.frame(chrom = x[1], pos = as.integer(x[2]), A2 = x[3], A1 = x[4])
    }))
    
    LD_data <- list(
        combined_LD_variants = ld_variant_id,
        combined_LD_matrix = LD,
        ref_panel = ref_panel
        )
    # check for duplicated positions
    variant_id = sumstats$variant_id

    alt_id <- variant_id |>
        as_tibble() |>
        separate(value, into = c("chrom", "pos", "A2", "A1"), sep = ":")|>
        mutate(alt_id = paste(chrom, pos, A1, A2, sep = ":")) |>
        pull(alt_id)


    error_list <- inner_join(as_tibble(variant_id), as_tibble(alt_id), by = join_by(value)) |>
        pull(value)
    
    # allele flip
    preprocess_results <- sumstats |>
        as_tibble() |>
        filter(!(variant_id %in% error_list)) |>
        rss_basic_qc(LD_data)
    sumstats <- preprocess_results$sumstats
  
    # remove outliers 
    qc_results <- summary_stats_qc(sumstats, LD_data, n = median(sumstats$n), 
                                   method = "rss_qc")
    
    return(qc_results)
}
#' Match Summary Statistics to LD matrix
#'
#' This function takes summary statistics for an entire chromosome and returns only the rows 
#' that fall within the position range of a given LD block.
#'
#' @param sumstats A dataframe of summary statistics including a variable named pos for the 
#'        position of each variant.
#' @param LD An LD matrix with column names in chrom:pos:A2:A1 format.
#'
#' @return A dataframe that is a subset of the provided summary statistics, with positions 
#'         aligning to the region of the LD matrix.
#'
#' @import stringr dplyr
#' @export
match_ld <- function(sumstats, LD){
    
    # get start and end position of LD region
    reg_start = head(colnames(LD), 1) |>
        str_extract("(?<=\\d:)\\d*(?=:)") |>
        as.numeric()
    reg_end = tail(colnames(LD), 1) |>
        str_extract("(?<=\\d:)\\d*(?=:)") |>
        as.numeric()
    
    # filter sumstats to within region
    sumstats_filtered <- sumstats |>
        filter(pos >= reg_start & pos <= reg_end)
    
    return(sumstats_filtered)
}
#' Format RSS 
#'
#' This function modifies summary statistics from Sun et al., selecting only necessary columns and naming 
#' the columns as needed by subsequent functions.
#'
#' @import dplyr
#' @export
rss_format <- function(df) {
    df = df |>
        rename(
            chrom = CHROM,
            pos = GENPOS,
            A2 = ALLELE0,
            A1 = ALLELE1,
            maf = A1FREQ,
            n = N,
            beta = BETA,
            se = SE,
            log10p = LOG10P) |>
        mutate(z = beta/se,
               variant_id = paste(chrom, pos, A2, A1, sep=":")
              ) |>
        select(-c(INFO, ID, TEST, CHISQ, EXTRA))
               
    return(df)
}
#' Untar and Extract Summary Statistics
#'
#' Given a path to a tarball containing summary statistics in separate .gz files, this function 
#' extracts the tarball and returns a list of dataframes, one for each chromosome.
#'
#' @param file_path The file path to the tarball of interest.
#' 
#' @return A list of summary statistic dataframes.
#'
#' @import stringr
#' @importFrom data.table fread
#' @export
untar_rss <- function(file_path){
    
    # Extract tarball into temporary directory
    temp_dir <- tempdir()
    untar(file_path, exdir = temp_dir)
    
    # Extract basename of tarball
    tarball <- basename(file_path)

    # List files in temporary diretory
    gz_files <- list.files(paste0(temp_dir, "/", str_extract(tarball, "^.*(?=\\.tar)"), "/", 
                                  sep=""), pattern = "\\.gz$", full.names = FALSE)

    # Initialize a list to store the matrices
    data_list <- list()

    # Loop through .gz files
    for (gz_file in gz_files) {
        chr_name <- basename(gz_file) |>
            str_extract("chr[\\d|X]*")
        
        data_list[[chr_name]] <- fread(
          paste0(temp_dir, "/", str_extract(tarball, "^.*(?=\\.tar)"), "/", gz_file, sep = "")) 
    }
    
    return(data_list)
}


#' Utility function for getting the statistical mode
#'
#' @export
stat_mode <- function(x) {
    unique_vals <- unique(x)
    unique_vals[which.max(tabulate(match(x, unique_vals)))]
}

#' Utility function for imputing NA values in a matrix
#'
#' @export
clean_matrix_na <- function(mat) {
  # Replace NA values in each column with the mode of that column
  mat_clean <- apply(mat, 2, function(col) {
    if (any(is.na(col))) {
      col[is.na(col)] <- stat_mode(col[!is.na(col)])
    }
    col
  })
  
  # Convert the result back to a matrix
  mat_clean <- as.matrix(mat_clean)
  
  # Return the cleaned matrix
  return(mat_clean)
}

#' Wrapper to obtain SuSiE weights
#'
#' This function imputes NA values into a genotype matrix based upon column modes, and then executes 
#' susie_weights
#'
#' @param y A vector of phenotypes
#' @param X A genotype matrix
#' @param susie_fit A susie fit object
#' 
#' @importFrom pecotmr susie_weights
#' @export
susie_weights_wrapper <- function(y, X, ...) {

    # Impute column mode to NA values
    # X_clean <- clean_matrix_na(X)
    print(sum(is.na(X)))
    # Call susie_weights function with the cleaned X
    weights <- X |>
        clean_matrix_na() |>
        susie_weights(y = y, X = _, ...)
    
    return(weights)
}

                               
#' Return weights from linear regression  
#' @export                               
lm_weights <- function(y, X){
    # Fit a simple linear regression
    fit = lm(y ~ X)
    
    # Return the coefficient for the variable
    return(fit$coefficients[2])
}