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


## Run spatialDE pipeline on any matrix-like object
#' @importFrom Matrix colSums
.spatialDE <- function(x, coordinates, verbose = FALSE) {
    sample_info <- data.frame(coordinates, total_counts = colSums(x))

    out <- basilisk::basiliskRun(
        env = spatialDE_env,
        fun = .run_spatialDE,
        x = x,
        sample_info = sample_info,
        verbose = verbose
    )
    out
}

.run_spatialDE <- function(x, sample_info, verbose = FALSE) {
    ## Normalization
    stabilized <- .naiveDE_stabilize(counts = x)
    regressed <- .naiveDE_regress_out(counts = stabilized, sample_info)

    coordinates <- sample_info[, c("x", "y")]
    .spatialDE_run(regressed, coordinates = coordinates, verbose = verbose)
}


#' @import methods
#' @export
#' @rdname spatialDE
setGeneric("spatialDE", function(x, ...) standardGeneric("spatialDE"))

#' @export
#' @rdname spatialDE
setMethod("spatialDE", "matrix", function(x, coordinates, ...) {
    .spatialDE(x = x, coordinates = coordinates, ...)
})

#' @export
#' @rdname spatialDE
#' @importFrom SummarizedExperiment assay
#' @importFrom SpatialExperiment spatialCoords spatialCoordsNames<-
setMethod("spatialDE", "SpatialExperiment",
    function(x, ..., assay_type = "counts") {
        ## Rename spatialCoords columns to "x", "y"
        spatialCoordsNames(x) <- c("x", "y")
        coordinates <- spatialCoords(x, as_df = TRUE)

        .spatialDE(x = assay(x, assay_type), coordinates = coordinates, ...)
    }
)
