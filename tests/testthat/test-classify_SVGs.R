set.seed(42)
mock <- mockSVG(10, 100, 10)
sample_info <- mock$coordinates
sample_info$total_counts <- colSums(mock$counts)

stabilized <- stabilize(mock$counts)
regressed <- regress_out(counts = stabilized, sample_info = sample_info)
de_results <- run(x = regressed, coordinates = mock$coordinates)

test_that("model_search() and spatial_patterns() return correct output", {
    out <- model_search(
        x = regressed, coordinates = mock$coordinates,
        de_results = de_results,
        qval_thresh = NULL   # disable filtering
    )
    expect_equal(nrow(out), nrow(de_results))
    expect_true(is.data.frame(out))
})

test_that("spatial_patterns() returns correct output", {
    n_patterns <- 1L
    sp <- spatial_patterns(
        x = regressed, coordinates = mock$coordinates,
        de_results = de_results, qval_thresh = NULL,  # disable filtering
        n_patterns = n_patterns, length = 1
    )
    pat_res <- sp$pattern_results
    pat <- sp$patterns
    expect_equal(nrow(pat_res), nrow(de_results))
    expect_equal(ncol(pat), n_patterns)
    expect_true(all(unique(pat_res$pattern) %in% colnames(pat)))
})

test_that("model_search() and spatial_patterns() break when necessary", {
    ## Check breaking errors (incompatible dimensions)
    expect_error(model_search(
        x = mock$counts, coordinates = mock$coordinates[1:3, ],
        de_results = de_results, qval_thresh = NULL
    ))

    ## Check for error when filtering everything out
    expect_error(model_search(
        x = mock$counts, coordinates = mock$coordinates,
        de_results = de_results, qval_thresh = 0
    ))
})

test_that("spatial_patterns() breaks when necessary", {
    ## Check breaking errors (incompatible dimensions)
    expect_error(spatial_patterns(
        x = mock$counts, coordinates = mock$coordinates[1:3, ],
        de_results = de_results, qval_thresh = NULL, n_patterns = 2L, length = 1
    ))

    ## Check for error when filtering everything out
    expect_error(spatial_patterns(
        x = mock$counts, coordinates = mock$coordinates,
        de_results = de_results, qval_thresh = 0,
        n_patterns = 2L, length = 1
    ))
})
