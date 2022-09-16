# TODO: add more details regarding model_search results in @return

#' Compare model fits with different models
#'
#' Classify DE genes to interpretable fitting classes.
#'
#' @inheritParams run
#' @param de_results `data.frame` resulting from [run()].
#' @param qval_thresh `numeric` scalar, specifying the q-value significance
#'   threshold to filter `de_results`. Only rows in `de_results` with
#'   `qval < qval_thresh` will be kept. To disable, set `qval_thresh = NULL`.
#'
#' @return `data.frame` of model_search results.
#'
#' @examples
#' set.seed(42)
#' mock <- mockSVG(size = 10, tot_genes = 300, de_genes = 10)
#' stabilized <- stabilize(mock$counts)
#' sample_info <- mock$coordinates
#' sample_info$total_counts <- colSums(mock$counts)
#' regressed <- regress_out(counts = stabilized, sample_info = sample_info)
#'
#' ## Run SpatialDE
#' de_results <- run(regressed, coordinates = mock$coordinates)
#'
#' ## Run model search
#' ms_results <- model_search(
#'     x = regressed,
#'     coordinates = mock$coordinates,
#'     de_results = de_results,
#'     qval_thresh = NULL
#' )
#'
#' @references Svensson, V., Teichmann, S. & Stegle, O. SpatialDE:
#'   identification of spatially variable genes. Nat Methods 15, 343–346 (2018).
#'   \url{https://doi.org/10.1038/nmeth.4636}
#'
#' @export
#' @importFrom checkmate assert_data_frame assert_names assert_matrix
#' @importFrom checkmate assert_flag
#' @importFrom basilisk basiliskRun basiliskStart
model_search <- function(x,
                         coordinates, de_results,
                         qval_thresh = 0.05,
                         verbose = FALSE) {
    assert_data_frame(coordinates, any.missing = FALSE)
    assert_names(colnames(coordinates), identical.to = c("x", "y"))
    assert_matrix(x, any.missing = FALSE)
    assert_data_frame(de_results, all.missing = FALSE)
    assert_names(colnames(de_results), must.include = "qval")
    assert_number(qval_thresh, null.ok = TRUE)
    assert_flag(verbose)

    ## Filter DE results
    if (!is.null(qval_thresh)) {
        de_results <- .filter_de_results(
            de_results = de_results, qval_thresh = qval_thresh
        )
    }
    
    proc <- basiliskStart(spatialDE_env, testload="scipy.optimize")
    .importPyModule(proc, !verbose, .set_fake_tqdm, .set_real_tqdm)

    .spatialDE_model_search(proc, x, coordinates, de_results)
    
    out <- basiliskRun(proc, function(store) {
        store$model_search
    }, persist=TRUE)
    
    out
}


#' @importFrom reticulate r_to_py
#' @importFrom basilisk basiliskRun
.spatialDE_model_search <- function(proc, x, coordinates, de_results) {
    basiliskRun(proc, function(x, coordinates, de_results, store) {
        spatialDE <- store$spatialDE
        X <- r_to_py(coordinates)
        
        ## Need to transpose counts for `spatialDE$model_search`
        res_py <- r_to_py(as.data.frame(t(x)))
        
        de_results_py <- r_to_py(de_results)
        
        out <- spatialDE$model_search(X, res_py, de_results_py)
        store$model_search <- out
        invisible(NULL)        
    }, x = x, coordinates = coordinates, de_results = de_results, persist=TRUE)
}



#' Group spatially variable genes into spatial patterns using automatic
#' expression histology (AEH)
#'
#' @inheritParams model_search
#' @param n_patterns `integer` The number of spatial patterns
#' @param length `numeric` The characteristic length scale of the clusters
#'
#' @return `list` of two dataframe (pattern_results, patterns):
#' `pattern_results` dataframe with pattern membership information for each
#' gene.
#' `patterns` the posterior mean underlying expression fro genes in given
#' spatial patterns.
#'
#' @examples
#' set.seed(42)
#' mock <- mockSVG(size = 10, tot_genes = 500, de_genes = 10)
#' stabilized <- stabilize(mock$counts)
#' sample_info <- mock$coordinates
#' sample_info$total_counts <- colSums(mock$counts)
#' regressed <- regress_out(counts = stabilized, sample_info = sample_info)
#'
#' ## Run SpatialDE
#' de_results <- run(x = regressed, coordinates = mock$coordinates)
#'
#' ## Run Spatial_patterns
#' sp <- spatial_patterns(
#'     x = regressed,
#'     coordinates = mock$coordinates,
#'     de_results = de_results,
#'     qval_thresh = NULL,
#'     n_patterns = 5, length = 1.5
#' )
#'
#' sp$pattern_results
#' sp$patterns
#'
#' @references
#' Svensson, V., Teichmann, S. & Stegle, O.
#' SpatialDE: identification of spatially variable genes.
#' Nat Methods 15, 343–346 (2018). \url{https://doi.org/10.1038/nmeth.4636}
#'
#' @export
#' @importFrom checkmate assert_data_frame assert_names assert_matrix
#' @importFrom checkmate assert_int assert_number assert_flag
#' @importFrom basilisk basiliskRun basiliskStart
spatial_patterns <- function(x, coordinates, de_results,
                             qval_thresh = 0.05,
                             n_patterns, length,
                             verbose = FALSE) {
    assert_data_frame(coordinates, any.missing = FALSE)
    assert_names(colnames(coordinates), identical.to = c("x", "y"))
    assert_matrix(x, any.missing = FALSE)
    assert_data_frame(de_results, all.missing = FALSE)
    assert_number(qval_thresh, null.ok = TRUE)
    assert_int(n_patterns, coerce = TRUE)
    assert_number(length)
    assert_flag(verbose)

    ## Filter de_results
    if (!is.null(qval_thresh)) {
        de_results <- .filter_de_results(
            de_results = de_results, qval_thresh = qval_thresh
        )
    }
    
    proc <- basiliskStart(spatialDE_env, testload="scipy.optimize")
    .importPyModule(proc, !verbose, .set_fake_tqdm, .set_real_tqdm)
    
    .spatialDE_spatial_patterns(proc, x, coordinates, de_results, n_patterns, 
                                length)

    out <- basiliskRun(proc, function(store) {
        store$spatial_pattern
    }, persist=TRUE)
    
    out
}


#' @importFrom reticulate r_to_py
#' @importFrom basilisk basiliskRun
.spatialDE_spatial_patterns <- function(proc, x, coordinates, de_results,
                                        n_patterns, length) {
    basiliskRun(proc, function(x, coordinates, de_results, n_patterns, length, 
                               store) {
        spatialDE <- store$spatialDE
        
        X <- r_to_py(coordinates)
        regr_out_py <- r_to_py(as.data.frame(t(x)))
        de_results_py <- r_to_py(de_results)
        
        spatterns <- spatialDE$spatial_patterns(
            X, regr_out_py, de_results_py,
            as.integer(n_patterns), length
        )
        
        out <- list(
            pattern_results = spatterns[[1]],
            patterns = spatterns[[2]]
        )
        
        store$spatial_pattern <- out        
        invisible(NULL)
    }, x=x, coordinates=coordinates, de_results = de_results,
    n_patterns = n_patterns, length = length, persist=TRUE)
}

## Helper to filter de_results and throw sensible error
.filter_de_results <- function(de_results, qval_thresh) {
    qval <- NULL
    if (qval_thresh < min(de_results$qval)) {
        stop(
            "Using `qval_thresh = ", qval_thresh, "` will filter out all genes",
            "\nConsider setting a less stringent threshold.",
            call. = FALSE
        )
    }
    subset(de_results, qval < qval_thresh)
}
