.spatialDE_dependencies <- c(
    "numpy==1.23.5",
    "patsy==0.5.3",
    "pandas==1.5.2"
)
#' @importFrom basilisk BasiliskEnvironment
spatialDE_env <- BasiliskEnvironment(
    envname = "env",
    pkgname = "spatialDE",
    packages = .spatialDE_dependencies,
    pip = c("SpatialDE==1.1.3", "NaiveDE==1.2.0")
)
