library(SpatialExperiment)

## Create mock SpatialExperiment object
spe <- mockSVG(size = 5, tot_genes = 100, de_genes = 10, return_SPE = TRUE)

counts <- counts(spe)
coordinates <- spatialCoords(spe)

test_that("spatialDE works with matrix input", {
    out <- spatialDE(counts, coordinates = coordinates, verbose = FALSE)
    expect_s3_class(out, "data.frame")
    expect_equal(nrow(out), nrow(spe))
})

test_that("spatialDE works with SpatialExperiment input", {
    ref <- spatialDE(counts, coordinates = coordinates, verbose = FALSE)
    out <- spatialDE(spe, assay_type = "counts", verbose = FALSE)
    ## Check if matching, ignore `time` column: this just measures computation
    ## time for a given gene
    out[["time"]] <- ref[["time"]] <- NULL
    expect_equal(out, ref, ignore_attr = TRUE)
})
