#' List available methods from SpatialDE
#'
#' Shows the available methods from the
#' [SpatialDE](https://github.com/Teichlab/SpatialDE) Python module.
#'
#' @export
#'
#' @examples
#' list_spatialDE_methods()
list_spatialDE_methods <- function() {
    names <- basilisk::basiliskRun(env = spatialDE_env, fun = function() {
        X <- reticulate::import("SpatialDE")
        names(X)
    })
    names
}

#' @note the function 'spatialDE$run' run with tqdm python library, which
#' allow to view a progress bar and it's not suitable for R.
#' TODO: hide progress bar.
#' 
#' @importFrom reticulate import
spatialDE_run <- function(coordinates, res) {
    spatialDE <- import("SpatialDE")
    
    X <- r_to_py(coordinates)
    res_py <- r_to_py(res)
    
    results <- spatialDE$run(X, res_py)
    return(results)
}

#' @importFrom reticulate import
spatialDE_model_search <- function(coordinates, res, de_results) {
    spatialDE <- import("SpatialDE")
    
    X <- r_to_py(coordinates)
    res_py <- r_to_py(res)
    de_results_py <- r_to_py(de_results)
    
    ms_results <- spatialDE$model_search(X, res_py, de_results_py)
    return(ms_results)
}
