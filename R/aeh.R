#' Automatic expression histology in **SpatialDE**
#'
#' Group spatially variable genes into spatial patterns using Automatic
#' Expression Histology, using the
#' [**SpatialDE**](https://github.com/Teichlab/SpatialDE) Python package.
#'
#' @param x A numeric `matrix` of counts where genes are rows and cells are columns.
#'
#'    Alternatively, a \linkS4class{SpatialExperiment} object.
#'
#' @param ... For the generic, arguments to pass to specific methods.
#' @param de_results `data.frame` resulting from [run()] or [spatialDE()].
#' @param coordinates A `data.frame` with sample coordinates. Each row is a
#'   sample, the columns with coordinates should be named 'x' and 'y'.
#'
#'   For the *SpatialExperiment* method, coordinates are taken from
#'   `spatialCoords(x)`.
#'
#' @param qval_thresh `numeric` scalar, specifying the q-value significance
#'   threshold to filter `de_results`. Only rows in `de_results` with
#'   `qval < qval_thresh` will be kept. To disable, set `qval_thresh = NULL`.
#' @param n_patterns `integer` The number of spatial patterns
#' @param length `numeric` The characteristic length scale of the clusters
#' @param assay_type A `character` string specifying the assay from `x` to use
#'   as input. Defaults to `"counts"`.
#' @param verbose A `logical` controlling the display of a progress bar from the
#'   Python package.
#'
#' @return A `list` of two `data.frame`s (pattern_results, patterns):
#'
#' * `pattern_results`: `data.frame` with pattern membership information for each
#' gene.
#'
#' * `patterns` the posterior mean underlying expression from genes in given
#' spatial patterns.
#'
#' @examples
#' ## Mock up a SpatialExperiment object wit 100 cells, 200 genes
#' set.seed(42)
#' spe <- mockSVG(size = 10, tot_genes = 200, de_genes = 20, return_SPE = TRUE)
#'
#' ## Run spatialDE
#' de_results <- spatialDE(spe)
#'
#' spatial_patterns <- spatialPatterns(spe, de_results = de_results,
#'     qval_thresh = NULL, n_patterns = 4L, length = 1.5,
#'     verbose = FALSE
#' )
#'
#' head(spatial_patterns$pattern_results)
#' head(spatial_patterns$patterns)
#'
#' @seealso
#' The individual steps performed by this function: [stabilize()],
#' [regress_out()] and [spatial_patterns()].
#'
#' @references
#' Svensson, V., Teichmann, S. & Stegle, O. SpatialDE: identification of
#' spatially variable genes. Nat Methods 15, 343–346 (2018).
#' \url{https://doi.org/10.1038/nmeth.4636}
#'
#' [**SpatialDE 1.1.3**](https://pypi.org/project/SpatialDE/1.1.3/): the version
#' of the Python package used under the hood.
#'
#' @author Davide Corso, Milan Malfait, Lambda Moses
#' @name spatialPatterns
NULL

#' @importFrom Matrix colSums
#' @importFrom checkmate assert_data_frame assert_number assert_int assert_flag
.spatialPatterns <- function(x, de_results, coordinates,
                             qval_thresh = 0.05, n_patterns, length,
                             verbose = FALSE) {
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
    
    sample_info <- data.frame(coordinates, total_counts = colSums(x))
    .run_spatial_patterns(x, sample_info, de_results, n_patterns, length,
                          verbose = FALSE)
}

#' @importFrom basilisk basiliskStart basiliskRun
.run_spatial_patterns <- function(x, sample_info, de_results,
                                  n_patterns, length,
                                  verbose = FALSE) {
    
    proc <- basiliskStart(spatialDE_env, testload="scipy.optimize")
    
    # Normalization
    ## Stabilize
    .naiveDE_stabilize(proc, x)
    stabilized <- basiliskRun(proc, function(store) {
        as.matrix(store$stabilized)
    }, persist=TRUE)
    
    ## Regress_out
    .naiveDE_regress_out(proc, stabilized, sample_info)
    regressed <- basiliskRun(proc, function(store) {
        as.matrix(store$regressed)
    }, persist=TRUE)
    
    coordinates <- sample_info[, c("x", "y")]
    
    .importPyModule(proc, !verbose, .set_fake_tqdm, .set_real_tqdm)
    .spatialDE_spatial_patterns(
        proc = proc, 
        x = regressed, coordinates = coordinates,
        de_results = de_results, n_patterns = n_patterns,
        length = length
    )
    
    out <- basiliskRun(proc, function(store) {
        store$spatial_pattern
    }, persist=TRUE)
    
    out
}

#' @import methods
#' @export
#' @rdname spatialPatterns
setGeneric("spatialPatterns",
    function(x, de_results, ...) standardGeneric("spatialPatterns"),
    signature = "x"
)

#' @export
#' @rdname spatialPatterns
setMethod("spatialPatterns", "matrix", .spatialPatterns)

#' @export
#' @rdname spatialPatterns
#' @importFrom SummarizedExperiment assay
#' @importFrom SpatialExperiment spatialCoords spatialCoordsNames<-
setMethod("spatialPatterns", "SpatialExperiment",
    function(x, de_results, qval_thresh = 0.05,
             n_patterns, length,
             assay_type = "counts", verbose = FALSE) {

        ## Rename spatialCoords columns to "x", "y"
        spatialCoordsNames(x) <- c("x", "y")
        coordinates <- as.data.frame(spatialCoords(x))

        .spatialPatterns(
            x = assay(x, assay_type), de_results = de_results,
            coordinates = coordinates, qval_thresh = qval_thresh,
            n_patterns = n_patterns, length = length, verbose = verbose
        )
    }
)
