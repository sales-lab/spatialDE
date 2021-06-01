#' Classify SVGenes to interpretable fitting classes
#'
#' Compare model fits with different models. 
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
#' @param verbose A `logical` controlling the display of a progress bar from the
#'   Python package.
#'
#' @return `data.frame` of model_search results.
#'
#' @examples
#' ## Mock up a SpatialExperiment object wit 100 cells, 200 genes
#' set.seed(42)
#' spe <- mockSVG(size = 10, tot_genes = 200, de_genes = 20, return_SPE = TRUE)
#'
#' ## Run spatialDE with S4 integration
#' de_results <- spatialDE(spe)
#' 
#' ## Run model search with S4 integration
#' model_search <- modelSearch(spe, assay_type = "counts", 
#' de_results = de_results, qval_thresh = NULL, verbose = FALSE)
#' 
#' 
#' @seealso
#' The individual steps performed by this function: [stabilize()],
#' [regress_out()] and [model_search()].
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
#' @name model-search
NULL

#' @importFrom Matrix colSums
.modelSearch <- function(counts_spe, coordinates_spe, de_results, 
                         verbose = FALSE) {
  sample_info <- data.frame(coordinates_spe, total_counts = colSums(counts_spe))
  
  out <- basilisk::basiliskRun(
    env = spatialDE_env,
    fun = .run_model_search,
    counts_spe=counts_spe, 
    sample_info=sample_info, 
    de_results=de_results,
    verbose=verbose
  )
  out
}

.run_model_search <- function(counts_spe, sample_info, de_results,
                              verbose = FALSE) {
  ## Normalization
  stabilized <- .naiveDE_stabilize(counts = counts_spe)
  regressed <- .naiveDE_regress_out(counts = stabilized, sample_info)
  
  coordinates <- sample_info[, c("x", "y")]
  .spatialDE_model_search(x = regressed, coordinates = coordinates,
                          de_results = de_results, verbose = verbose)
}

#' 
#' @import methods
#' @export
#' @rdname model-search
setGeneric("modelSearch", function(x, ...) 
  standardGeneric("modelSearch"))

#' @export
#' @rdname model-search
#' @importFrom SummarizedExperiment assay
#' @importFrom SpatialExperiment spatialCoords spatialCoordsNames
#' @importFrom checkmate assert_data_frame assert_number assert_flag
setMethod("modelSearch", "SpatialExperiment",
  function(x, assay_type = "counts", de_results, qval_thresh=0.05, 
           verbose = FALSE) {
    assert_data_frame(de_results, all.missing = FALSE)
    assert_number(qval_thresh, null.ok = TRUE)
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
    
    .modelSearch(counts_spe=counts_spe, coordinates_spe=coordinates_spe,
                 de_results=de_results, verbose = FALSE)
  }
)
