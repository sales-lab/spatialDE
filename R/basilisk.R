#' @importFrom basilisk BasiliskEnvironment
spatialDE_env <- BasiliskEnvironment(
    envname = "env", pkgname = "spatialDE",
    packages = character(0), pip = c("SpatialDE==1.1.3", "NaiveDE==1.2.0")
)
