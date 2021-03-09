#' Find spatially variable genes with **SpatialDE**
#'
#' Identify genes that significantly depend on spatial coordinates with the
#' [**SpatialDE**](https://github.com/Teichlab/SpatialDE) Python package.
#'
#' @param x A numeric `matrix` of counts where genes are rows and cells are columns.
#'
#'    Alternatively, a \linkS4class{SpatialExperiment} object.
#'
#' @param ... For the generic, arguments to pass to specific methods.
#' @param coordinates A `data.frame` with sample coordinates. Each row is a
#'   sample, the columns with coordinates should be named 'x' and 'y'.
#'
#'   For the *SpatialExperiment* method, coordinates are taken from
#'   `spatialCoords(x)`.
#'
#' @param assay_type A `character` string specifying the assay from `x` to use
#'   as input. Defaults to `"counts"`.
#' @param verbose A `logical` controlling the display of a progress bar from the
#'   Python package.
#'
#' @return A `data.frame` with DE results where each row is a gene and columns
#'   contain relevant statistics.
#'
#'   The most important columns are:
#'
#'   * `g`: the name of the gene
#'   * `pval`: the p-value for spatial differential expression
#'   * `qval`: the q-value, indicating significance after correcting for
#'   multiple testing
#'   * `l`: A parameter indicating the distance scale a gene changes expression
#'   over
#'
#' @examples
#' ## Mock up a SpatialExperiment object wit 100 cells, 200 genes
#' set.seed(42)
#' spe <- mockSVG(size = 10, tot_genes = 200, de_genes = 20, return_SPE = TRUE)
#'
#' ## Run spatialDE
#' de_results <- spatialDE(spe)
#' head(de_results)
#'
#' @seealso
#' The individual steps performed by this function: [stabilize()],
#' [regress_out()] and [run()].
#'
#' For further analysis of the DE results:
#' [model_search()] and [spatial_patterns()].
#'
#' @references
#' Svensson, V., Teichmann, S. & Stegle, O. SpatialDE: identification of
#' spatially variable genes. Nat Methods 15, 343–346 (2018).
#' \url{https://doi.org/10.1038/nmeth.4636}
#'
#' [**SpatialDE 1.1.3**](https://pypi.org/project/SpatialDE/1.1.3/): the version
#' of the Python package used under the hood.
#'
#' @author Davide Corso, Milan Malfait
#' @name spatialDE
NULL

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
setMethod("spatialDE", "matrix", .spatialDE)

#' @export
#' @rdname spatialDE
#' @importFrom SummarizedExperiment assay
#' @importFrom SpatialExperiment spatialCoords spatialCoordsNames<-
setMethod("spatialDE", "SpatialExperiment",
    function(x, assay_type = "counts", verbose = FALSE) {
        ## Rename spatialCoords columns to "x", "y"
        spatialCoordsNames(x) <- c("x", "y")
        coordinates <- spatialCoords(x, as_df = TRUE)

        .spatialDE(
            x = assay(x, assay_type),
            coordinates = coordinates,
            verbose = verbose
        )
    }
)
