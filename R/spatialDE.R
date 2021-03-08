# TODO: add more details regarding DE results format in @return

#' Perform SpatialDE test
#'
#' Wraps the `run` function from the
#' [SpatialDE](https://github.com/Teichlab/SpatialDE) Python package.
#'
#' @param x `matrix`-like object of normalized counts. E.g. resulting from
#'          [regress_out()].
#' @param coordinates `data.frame` with sample coordinates.
#' Each row is a sample, the columns with coordinates must be named 'x' and 'y'.
#' @param verbose `logical` controlling the display of the progress bar.
#'
#' @return `data.frame` with DE results.
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
#' de_results <- run(regressed, coordinates = mock$coordinates)
#'
#' @references
#' Svensson, V., Teichmann, S. & Stegle, O.
#' SpatialDE: identification of spatially variable genes.
#' Nat Methods 15, 343–346 (2018). \url{https://doi.org/10.1038/nmeth.4636}
#'
#' @export
#' @importFrom checkmate assert_data_frame assert_names assert_matrix
#' @importFrom checkmate assert_flag
run <- function(x, coordinates, verbose = FALSE) {
    assert_data_frame(coordinates, any.missing = FALSE)
    assert_names(colnames(coordinates), identical.to = c("x", "y"))
    assert_matrix(x, any.missing = FALSE)
    assert_flag(verbose)

    out <- basilisk::basiliskRun(
        env = spatialDE_env,
        fun = .spatialDE_run,
        x = x,
        coordinates = coordinates,
        verbose = verbose
    )
    out
}

#' @importFrom reticulate r_to_py
.spatialDE_run <- function(x, coordinates, verbose) {
    spatialDE <- .importPyModule(!verbose)

    X <- r_to_py(coordinates)
    ## Need to transpose counts for `spatialDE$run`
    res_py <- r_to_py(as.data.frame(t(x)))

    spatialDE$run(X, res_py)
}


# TODO: add more details regarding model_search results in @return

#' Compare model fits with different models
#'
#' Classify DE genes to interpretable fitting classes.
#'
#' @inheritParams run
#' @param de_results `data.frame` resulting from [run()].
#' @param filter `logical`, whether to filter `de_results` based on
#'   `qval_thresh`. Default: `TRUE`
#' @param qval_thresh `numeric` scalar, specifying the q-value significance
#'   threshold to filter `de_results`. Only rows in `de_results` with
#'   `qval < qval_thresh` will be kept. Ignored if `filter = FALSE`.
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
#'     filter = FALSE
#' )
#'
#' @references Svensson, V., Teichmann, S. & Stegle, O. SpatialDE:
#'   identification of spatially variable genes. Nat Methods 15, 343–346 (2018).
#'   \url{https://doi.org/10.1038/nmeth.4636}
#'
#' @export
#' @importFrom checkmate assert_data_frame assert_names assert_matrix
#' @importFrom checkmate assert_flag
model_search <- function(x,
                         coordinates, de_results,
                         filter = TRUE,
                         qval_thresh = 0.05,
                         verbose = FALSE) {
    assert_data_frame(coordinates, any.missing = FALSE)
    assert_names(colnames(coordinates), identical.to = c("x", "y"))
    assert_matrix(x, any.missing = FALSE)
    assert_data_frame(de_results, all.missing = FALSE)
    assert_names(colnames(de_results), must.include = "qval")
    assert_flag(verbose)

    ## Filter DE results
    if (filter) {
        de_results <- .filter_de_results(
            de_results = de_results, qval_thresh = qval_thresh
        )
    }

    out <- basilisk::basiliskRun(
        env = spatialDE_env,
        fun = .spatialDE_model_search,
        x = x,
        coordinates = coordinates,
        de_results = de_results, verbose = verbose
    )
    out
}


#' @importFrom reticulate r_to_py
.spatialDE_model_search <- function(x, coordinates, de_results,
                                    verbose) {
    spatialDE <- .importPyModule(!verbose)

    X <- r_to_py(coordinates)

    ## Need to transpose counts for `spatialDE$model_search`
    res_py <- r_to_py(as.data.frame(t(x)))

    de_results_py <- r_to_py(de_results)

    spatialDE$model_search(X, res_py, de_results_py)
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
#'     filter = FALSE,
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
spatial_patterns <- function(x, coordinates, de_results,
                             filter = TRUE,
                             qval_thresh = 0.05,
                             n_patterns, length,
                             verbose = FALSE) {
    assert_data_frame(coordinates, any.missing = FALSE)
    assert_names(colnames(coordinates), identical.to = c("x", "y"))
    assert_matrix(x, any.missing = FALSE)
    assert_data_frame(de_results, all.missing = FALSE)
    assert_int(n_patterns, coerce = TRUE)
    assert_number(length)
    assert_flag(verbose)

    ## Filter de_results
    if (filter) {
        de_results <- .filter_de_results(
            de_results = de_results, qval_thresh = qval_thresh
        )
    }

    out <- basilisk::basiliskRun(
        env = spatialDE_env,
        fun = .spatialDE_spatial_patterns,
        x = x,
        coordinates = coordinates,
        de_results = de_results,
        n_patterns = n_patterns,
        length = length,
        verbose = verbose
    )
    out
}


#' @importFrom reticulate r_to_py
.spatialDE_spatial_patterns <- function(x, coordinates, de_results,
                                        n_patterns, length,
                                        verbose) {
    spatialDE <- .importPyModule(!verbose)

    X <- r_to_py(coordinates)
    regr_out_py <- r_to_py(as.data.frame(t(x)))
    de_results_py <- r_to_py(de_results)

    spatterns <- spatialDE$spatial_patterns(
        X, regr_out_py, de_results_py,
        as.integer(n_patterns), length
    )
    list(
        pattern_results = spatterns[[1]],
        patterns = spatterns[[2]]
    )
}

## Helper to filter de_results and throw sensible error
.filter_de_results <- function(de_results, qval_thresh) {
    if (qval_thresh < min(de_results$qval)) {
        stop(
            "Using `qval_thresh = ", qval_thresh, "` will filter out all genes",
            "\nConsider setting a less stringent threshold.",
            call. = FALSE
        )
    }
    subset(de_results, qval < qval_thresh)
}
