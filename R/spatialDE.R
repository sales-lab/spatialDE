#' Wrapper for SpatialDE.run python function
#'
#' @param coordinates `data.frame` with sample coordinates.
#' Each row is a sample, the columns with coordinates must be named 'x' and 'y'.
#'
#' @param regressed_counts `data.frame` or `matrix` resulting from
#' `naiveDE_regress_out`
#'
#' @note both functions of SpatialDE run with tqdm python library, which
#' allow to view a progress bar and it's not suitable for R.
#' TODO: hide prints of progress bar
#'
#' @export
run_spatialDE <- function(coordinates, regressed_counts) {
    out <- basilisk::basiliskRun(
        env = spatialDE_env,
        fun = .spatialDE_run,
        coordinates = coordinates, regressed_counts = regressed_counts
    )
    out
}


#' @importFrom reticulate import r_to_py
.spatialDE_run <- function(coordinates, regressed_counts) {
    spatialDE <- import("SpatialDE")

    X <- r_to_py(coordinates)
    ## Need to transpose counts for `spatialDE$run`
    res_py <- r_to_py(as.data.frame(t(regressed_counts)))

    spatialDE$run(X, res_py)
}



#' Wrapper for SpatialDE.model_search python function
#'
#' @param coordinates `data.frame` with sample coordinates.
#' Each row is a sample, the columns with coordinates must be named 'x' and 'y'.
#'
#' @param regressed_counts `data.frame` resulting from `naiveDE_regress_out`
#'
#' @param de_results `data.frame` resulting from `spatialDE_run` filtered based
#' on qvalue < threshold (e.g. qvalue < 0.05)
#'
#' @export
run_model_search <- function(coordinates, regressed_counts, de_results) {
    out <- basilisk::basiliskRun(
        env = spatialDE_env,
        fun = .spatialDE_model_search,
        coordinates = coordinates, regressed_counts = regressed_counts,
        de_results = de_results
    )
    out
}


#' @importFrom reticulate import r_to_py
.spatialDE_model_search <- function(coordinates, regressed_counts, de_results) {
    spatialDE <- import("SpatialDE")

    X <- r_to_py(coordinates)

    ## Need to transpose counts for `spatialDE$model_search`
    res_py <- r_to_py(as.data.frame(t(regressed_counts)))

    de_results_py <- r_to_py(de_results)

    spatialDE$model_search(X, res_py, de_results_py)
}
