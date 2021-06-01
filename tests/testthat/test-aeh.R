library(SpatialExperiment)

## Create mock SpatialExperiment object
spe <- mockSVG(size = 5, tot_genes = 100, de_genes = 10, return_SPE = TRUE)

counts <- counts(spe)
coordinates <- spatialCoords(spe)

de_results <- spatialDE(spe, assay_type = "counts", verbose = FALSE)

test_that("spatialPatterns works with matrix input", {
    out <- spatialPatterns(counts,
        de_results = de_results, coordinates = coordinates,
        qval_thresh = NULL, n_patterns = 4L, length = 1.5
    )
    expect_type(out, "list")
    expect_s3_class(out$patterns, "data.frame")
    expect_s3_class(out$pattern_results, "data.frame")

    expect_equal(nrow(out$pattern_results), nrow(spe))
    expect_equal(nrow(out$patterns), ncol(spe))
})

test_that("spatialPatterns works with SpatialExperiment input", {
    # ref <- spatialPatterns(counts,
    #     de_results = de_results, coordinates = coordinates,
    #     qval_thresh = NULL, n_patterns = 4L, length = 1.5
    # )
    out <- spatialPatterns(spe, de_results = de_results,
        qval_thresh = NULL, n_patterns = 4L, length = 1.5
    )
    expect_type(out, "list")
    expect_s3_class(out$patterns, "data.frame")
    expect_s3_class(out$pattern_results, "data.frame")

    expect_equal(nrow(out$pattern_results), nrow(spe))
    expect_equal(nrow(out$patterns), ncol(spe))

    ## Note that there is some stochasticity in the spatialPatterns results so
    ## we cannot directly compare the results to a reference
})
