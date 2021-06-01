#' Automatic expression histology in **SpatialDE**
#'
#' Group spatially variable genes into spatial patterns using Automatic 
#' Expression Histology. 
#' [**SpatialDE**](https://github.com/Teichlab/SpatialDE) Python package.
#'
#' @param x \linkS4class{SpatialExperiment} object.
#'
#' @param ... For the generic, arguments to pass to specific methods.
#' @param assay_type A `character` string specifying the assay from `x` to use
#'   as input. Defaults to `"counts"`.
#' @param de_results `data.frame` resulting from [run()] or [spatialDE()].
#' @param qval_thresh `numeric` scalar, specifying the q-value significance
#'   threshold to filter `de_results`. Only rows in `de_results` with
#'   `qval < qval_thresh` will be kept. To disable, set `qval_thresh = NULL`.
#' @param n_patterns `integer` The number of spatial patterns
#' @param length `numeric` The characteristic length scale of the clusters
#' @param verbose A `logical` controlling the display of a progress bar from the
#'   Python package.
#'
#' @return `list` of two dataframe (pattern_results, patterns):
#' `pattern_results` dataframe with pattern membership information for each
#' gene.
#' `patterns` the posterior mean underlying expression fro genes in given
#' spatial patterns.
#'
#' @examples
#' ## Mock up a SpatialExperiment object wit 100 cells, 200 genes
#' set.seed(42)
#' spe <- mockSVG(size = 10, tot_genes = 200, de_genes = 20, return_SPE = TRUE)
#'
#' ## Run spatialDE with S4 integration
#' de_results <- spatialDE(spe)
#' 
#' spatial_patterns <- spatialPatterns(spe, assay_type = "counts", 
#' de_results = de_results, qval_thresh = NULL, n_patterns = 4L, length = 1.5, 
#' verbose = FALSE)
#' 
#' head(de_results)
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
#' @author Davide Corso, Milan Malfait
#' @name spatial-patterns
NULL

#' @importFrom Matrix colSums
.spatialPatterns <- function(counts_spe, coordinates_spe, de_results, 
                             n_patterns, length, verbose = FALSE) {
  sample_info <- data.frame(coordinates_spe, total_counts = colSums(counts_spe))
  
  out <- basilisk::basiliskRun(
    env = spatialDE_env,
    fun = .run_spatial_patterns,
    counts_spe=counts_spe, 
    sample_info=sample_info, 
    de_results=de_results, 
    n_patterns=n_patterns, 
    length=length, 
    verbose=verbose
  )
  out
}

.run_spatial_patterns <- function(counts_spe, sample_info, 
                                  de_results, n_patterns, length, 
                                  verbose = FALSE) {
  ## Normalization
  stabilized <- .naiveDE_stabilize(counts = counts_spe)
  regressed <- .naiveDE_regress_out(counts = stabilized, sample_info)

  coordinates <- sample_info[, c("x", "y")]
  .spatialDE_spatial_patterns(x=regressed, coordinates=coordinates, 
                              de_results=de_results, n_patterns=n_patterns, 
                              length=length, verbose=verbose)
}

#' 
#' @import methods
#' @export
#' @rdname spatial-patterns
setGeneric("spatialPatterns", function(x, ...) 
  standardGeneric("spatialPatterns"))

#' @export
#' @rdname spatial-patterns
#' @importFrom SummarizedExperiment assay
#' @importFrom checkmate assert_data_frame assert_number assert_int assert_flag
#' @importFrom SpatialExperiment spatialCoords spatialCoordsNames
setMethod("spatialPatterns", "SpatialExperiment",
  function(x, assay_type = "counts", de_results, qval_thresh=0.05, 
            n_patterns, length, verbose = FALSE) {
    assert_data_frame(de_results, all.missing = FALSE)
    assert_number(qval_thresh, null.ok = TRUE)
    assert_int(n_patterns, coerce = TRUE)
    assert_number(length)
    assert_flag(verbose)
    
    ## Rename spatialCoords columns to "x", "y"
    spatialCoordsNames(x) <- c("x", "y")
    coordinates_spe <- as.data.frame(spatialCoords(x))
    counts_spe <- assay(x, assay_type)
    
    ## Filter de_results
    if (!is.null(qval_thresh)) {
      de_results <- .filter_de_results(
        de_results = de_results, qval_thresh = qval_thresh
      )
    }
    
    .spatialPatterns(counts_spe=counts_spe, coordinates_spe=coordinates_spe, 
                     de_results=de_results, n_patterns=n_patterns, 
                     length=length, verbose = FALSE)
  }
)
