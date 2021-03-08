

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
