#' Run multiple pQTL weight methods
#'
#' Applies specified methods to a list of genetic and phenotypic information. Works on either
#' regression summary statistics (RSS) or on individual-level data. This function (if directed)
#' utilizes parallel processing to handle multiple methods.
#'
#' @param y A vector of protein expression levels.
#' @param X A matrix of genotypes, where each row represents an individual and each column
#'        represents a SNP.
#' @param sumstats A dataframe with marker summary statistics. Required: beta coefficient (beta),
#'        standard error of the beta coefficient (se), GWAS sample size (n). Optional: rsids,
#'        alleles (A1 and A2), minor allele frequency (maf).
#' @param LD An LD matrix corresponding to the SNPs provided in the stat dataframe.
#' @param weight_methods A list of methods and their specific arguments, formatted as
#'        list(method1 = method1_args, method2 = method2_args). These methods will be applied
#'        to the data.
#' @param num_threads The number of threads to use for parallel processing.
#'        If set to -1, the function uses all available cores.
#'        If set to 0 or 1, no parallel processing is performed.
#'        If set to 2 or more, parallel processing is enabled with that many threads.
#' @param seed An optional seed for the reproducibility of the analysis.
#'
#' @return A list where each element is named after a method and contains the weight matrix
#'         produced by that method.
#'
#' @importFrom future availableCores
#' @export
pqtl_weights <- function(y = NULL, X = NULL, sumstats = NULL, LD = NULL,
                         weight_methods, num_threads = 1, seed = NULL) {
    # define approach
    if (!(is.null(y) | is.null(X))) {
        approach <- "individual"
    } else if (!(is.null(sumstats) | is.null(LD))) {
        approach <- "rss"
    } else {
        stop("Must provide either individual data or summary statistics.")
    }

    # determine number of cores to use
    num_cores <- ifelse(num_threads == -1, availableCores(), num_threads)
    num_cores <- min(num_cores, availableCores())

    # get weights
    if (approach == "individual") {
        weights_list <- compute_pqtl_weights(
            y = y, X = as.matrix(X),
            weight_methods = weight_methods,
            num_cores = num_cores, seed = seed
        )
    } else if (approach == "rss") {
        weights_list <- compute_pqtl_rss_weights(
            sumstats = sumstats, LD = LD, weight_methods = weight_methods,
            num_cores = num_cores, seed = seed
        )
    }

    # combine into a single matrix
    results <- do.call(cbind, weights_list)
    colnames(results) <- names(weight_methods)

    return(results)
}


#' Compute pQTL weights on individual data - WIP
#'
#' This function provides pQTL weights for SNPs given individual data and a set of methods
#'
#' @param y A vector of protein expression levels.
#' @param X A matrix of genotypes, where each row represents an individual and each column represents a SNP
#' @param weight_methods A list of methods and their specific arguments, formatted as
#'        list(method1 = method1_args, method2 = method2_args). These methods will be applied to the data.
#' @param num_cores The number of cores to use for parallel processing.
#' @param seed An optional seed for the reproducibility of the analysis.
#'
#' @return A list where each element is named after a method and contains the weight matrix produced by that method.
#'
#' @importFrom future plan multisession
#' @importFrom furrr future_map furrr_options
#' @importFrom purrr map
#' @export
compute_pqtl_weights <- function(y, X, weight_methods, num_cores, seed = NULL) {
    # define function
    compute_method_weights <- function(method_name) {
        args <- weight_methods[[method_name]]

        X_imputed <- apply(X, 2, function(col) {
            if (any(is.na(col))) {
                col[is.na(col)] <- stat_mode(col[!is.na(col)])
            }
            return(col)
        })
        valid_columns <- apply(
            X_imputed, 2,
            function(col) sd(col, na.rm = TRUE) != 0
        )

        X_filtered <- X_imputed[, valid_columns, drop = FALSE]
        
        # Initialize y with zeros to avoid NA
        if (!is.null(seed)) set.seed(seed)

        weights_matrix <- matrix(0, nrow = ncol(X_filtered), ncol = 1)
        weights_vector <- do.call(method_name, c(list(X = X_filtered, y = y), args))
        
        # Name and return weights
        if (method_name == "lm_weights") {
            names(weights_vector) <- colnames(X_filtered)
        }
        weights_matrix[, 1] <- weights_vector
        rownames(weights_matrix) <- names(weights_vector)
        
        return(weights_matrix)
    }

    # Execute function with or without multiprocessing
    if (num_cores >= 2) {
        plan(multisession, workers = num_cores)
        weights_list <- names(weight_methods) |> future_map(compute_method_weights,
            .options = furrr_options(seed = seed)
        )
    } else {
        weights_list <- names(weight_methods) |> map(compute_method_weights)
    }

    return(weights_list)
}


#' Compute pQTL weights on summary statistics
#'
#' This function provides pQTL weights for SNPs given summary statistics and a set of methods
#'
#' @param sumstats A dataframe with marker summary statistics. Required: beta coefficient (beta),
#'        standard error of the beta coefficient (se), GWAS sample size (n). Optional: rsids,
#'        alleles (A1 and A2), minor allele frequency (maf).
#' @param LD An LD matrix corresponding to the SNPs provided in the stat dataframe.
#' @param weight_methods A list of methods and their specific arguments, formatted as
#'        list(method1 = method1_args, method2 = method2_args). These methods will be applied to the data.
#' @param num_cores The number of cores to use for parallel processing.
#' @param seed An optional seed for the reproducibility of the analysis.
#'
#' @importFrom future plan multisession
#' @importFrom furrr future_map furrr_options
#' @importFrom purrr map
#' @export
compute_pqtl_rss_weights <- function(sumstats, LD, weight_methods, num_cores, seed = NULL) {
    # define function
    compute_method_weights <- function(method_name) {
        # Load pecotmr library with error handling
        tryCatch({
            library(pecotmr, quietly = TRUE)
        }, error = function(e) {
            stop("Failed to load pecotmr package: ", e$message)
        })
        args <- weight_methods[[method_name]]

        # Initialize Y with zeros to avoid NA
        if (!is.null(seed)) set.seed(seed)

        weights_matrix <- matrix(0, nrow = nrow(sumstats), ncol = 1)
        weights_vector <- do.call(method_name, c(list(sumstats = sumstats, LD = LD), args))
        weights_matrix[, 1] <- weights_vector

        return(weights_matrix)
    }

    # execute function with or without multiprocessing
    if (num_cores >= 2) {
        plan(multisession, workers = num_cores)
        weights_list <- names(weight_methods) |> future_map(compute_method_weights,
            .options = furrr_options(seed = seed)
        )
    } else {
        weights_list <- names(weight_methods) |> map(compute_method_weights)
    }

    return(weights_list)
}
