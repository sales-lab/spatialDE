.spatialDE_dependencies <- c(
    "click==7.1.2",
    "numpy==1.20.0",
    "pandas==1.2.1",
    "patsy==0.5.1",
    "pip==21.0.1",
    "python-dateutil==2.8.1",
    "scipy==1.6.0",
    "six==1.15.0",
    "tqdm==4.56.0",
    "wheel==0.36.2"
)
#' @importFrom basilisk BasiliskEnvironment
spatialDE_env <- BasiliskEnvironment(
    envname = "env", pkgname = "spatialDE",
    packages = .spatialDE_dependencies,
    pip = c("SpatialDE==1.1.3", "NaiveDE==1.2.0")
)
