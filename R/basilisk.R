.spatialDE_dependencies <- c(
    "python==3.11",
    "numpy==1.26.4",
    "scipy==1.11.4",
    "patsy==1.0.1",
    "pandas==1.5.3"
)
#' @importFrom basilisk BasiliskEnvironment
spatialDE_env <- BasiliskEnvironment(
    envname = "env",
    pkgname = "spatialDE",
    packages = .spatialDE_dependencies,
    pip = c("SpatialDE==1.1.3", "NaiveDE==1.2.0")
)
