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
