# TODO: add more details regarding DE results format in @return

#' Perform SpatialDE test
#'
#' Wraps the `run` function from the
#' [SpatialDE](https://github.com/Teichlab/SpatialDE) Python package.
#'
#' @param coordinates `data.frame` with sample coordinates.
#' Each row is a sample, the columns with coordinates must be named 'x' and 'y'.
#'
#' @param regressed_counts `data.frame` or `matrix` resulting from
#' [regress_out()].
#'
#' @param verbose `bool` controlling the display of the progress bar.
#'
#' @return `data.frame` with DE results.
#'
#' @examples
#' \dontrun{
#' ncells <- 100
#' ngenes <- 1000
#' counts <- matrix(rpois(ncells * ngenes, lambda = 3),
#'                  nrow = ngenes, ncol = ncells)
#'
#' ## Provide total counts for library size normalization and coordinates
#' sample_info <- data.frame(total_counts = colSums(counts))
#' coordinates <- data.frame(x = rnorm(ncells), y = rnorm(ncells))
#'
#' stabilized <- stabilize(counts)
#' regressed <- regress_out(sample_info, stabilized)
#'
#' ## Run SpatialDE
#' de_results <- run(coordinates, regressed)
#' }
#' @export
run <- function(coordinates, regressed_counts, verbose = FALSE) {
    out <- basilisk::basiliskRun(
        env = spatialDE_env,
        fun = .spatialDE_run,
        coordinates = coordinates, regressed_counts = regressed_counts,
        verbose = verbose
    )
    out
}

# NOTE: both functions of SpatialDE run with tqdm python library, which
# allow to view a progress bar and it's not suitable for R.
# TODO: hide prints of progress bar

#' @importFrom reticulate r_to_py
.spatialDE_run <- function(coordinates, regressed_counts, verbose) {
    spatialDE <- .importPyModule(!verbose)

    X <- r_to_py(coordinates)
    ## Need to transpose counts for `spatialDE$run`
    res_py <- r_to_py(as.data.frame(t(regressed_counts)))

    spatialDE$run(X, res_py)
}


# TODO: add more details regarding model_search results in @return

#' Compare model fits with different models
#'
#' Classify DE genes to interpretable fitting classes.
#'
#' @param coordinates `data.frame` with sample coordinates.
#' Each row is a sample, the columns with coordinates must be named 'x' and 'y'.
#'
#' @param regressed_counts `data.frame` resulting from [regress_out()]
#'
#' @param de_results `data.frame` resulting from [run()] filtered
#' based on `qvalue < threshold` (e.g. `qvalue < 0.05`)
#'
#' @param verbose `bool` controlling the display of the progress bar.
#'
#' @return `data.frame` of model_search results.
#'
#' @examples
#' \dontrun{
#' ncells <- 100
#' ngenes <- 1000
#' counts <- matrix(rpois(ncells * ngenes, lambda = 3),
#'                  nrow = ngenes, ncol = ncells)
#'
#' ## Provide total counts for library size normalization and coordinates
#' sample_info <- data.frame(total_counts = colSums(counts))
#' coordinates <- data.frame(x = rnorm(ncells), y = rnorm(ncells))
#'
#' stabilized <- stabilize(counts)
#' regressed <- regress_out(sample_info, stabilized)
#'
#' ## Run SpatialDE
#' de_results <- run(coordinates, regressed)
#'
#' ## Run model search
#' ms_results <- model_search(coordinates, regressed, de_results)
#' }
#' @export
model_search <- function(coordinates, regressed_counts, de_results,
                             verbose = FALSE) {
    out <- basilisk::basiliskRun(
        env = spatialDE_env,
        fun = .spatialDE_model_search,
        coordinates = coordinates, regressed_counts = regressed_counts,
        de_results = de_results, verbose = verbose
    )
    out
}


#' @importFrom reticulate r_to_py
.spatialDE_model_search <- function(coordinates, regressed_counts, de_results,
                                    verbose) {
    spatialDE <- .importPyModule(!verbose)

    X <- r_to_py(coordinates)

    ## Need to transpose counts for `spatialDE$model_search`
    res_py <- r_to_py(as.data.frame(t(regressed_counts)))

    de_results_py <- r_to_py(de_results)

    spatialDE$model_search(X, res_py, de_results_py)
}
