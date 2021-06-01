library(SpatialExperiment)

## Create mock SpatialExperiment object
spe <- mockSVG(size = 5, tot_genes = 100, de_genes = 10, return_SPE = TRUE)

counts <- counts(spe)
coordinates <- spatialCoords(spe)

de_results <- spatialDE(spe, assay_type = "counts", verbose = FALSE)

test_that("modelSearch works with matrix input", {
    out <- modelSearch(counts,
        de_results = de_results,
        coordinates = coordinates, qval_thresh = NULL
    )
    expect_s3_class(out, "data.frame")
    expect_equal(nrow(out), nrow(spe))
})

test_that("modelSearch works with SpatialExperiment input", {
    ref <- modelSearch(counts,
        de_results = de_results,
        coordinates = coordinates, qval_thresh = NULL
    )
    out <- modelSearch(spe,
        assay_type = "counts",
        de_results = de_results, qval_thresh = NULL,
        verbose = FALSE
    )
    ## Check if matching, ignore `time` column: this just measures computation
    ## time for a given gene
    out[["time"]] <- ref[["time"]] <- NULL
    expect_equal(out, ref, ignore_attr = TRUE)
})
