.spatialDE_dependencies <- c(
    #"click==8.0.1",
    #"numpy==1.21.0",
    "python==3.8.5",
    "pandas==0.25.3",
    "patsy==0.5.1",
    "pip==21.0.1"#,
    #"python-dateutil==2.8.1",
    #"scipy==1.7.0",
    #"six==1.16.0",
    #"tqdm==4.61.1",
    #"wheel==0.36.2"
)
#' @importFrom basilisk BasiliskEnvironment
spatialDE_env <- BasiliskEnvironment(
    envname = "env", pkgname = "spatialDE",
    packages = .spatialDE_dependencies,
    pip = c("SpatialDE==1.1.3", "NaiveDE==1.2.0")
)
