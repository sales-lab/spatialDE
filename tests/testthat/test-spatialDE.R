ncells <- 10
ngenes <- 100

## Note: `stabilize()` returns all NA values for some count inputs
## Setting a seed to get the same counts for which stabilize() works for now
set.seed(1)
counts <- matrix(stats::rpois(ncells * ngenes, lambda = 3),
                 nrow = ngenes, ncol = ncells)
sample_info <- data.frame(x = stats::rnorm(ncells), y = stats::rnorm(ncells),
                          total_counts = colSums(counts))
coordinates <- sample_info[, c("x", "y")]

test_that("Wrapper functions work", {
    stabilized <- stabilize(counts = counts)
    expect_false(all(is.na(stabilized)))

    regressed <- regress_out(counts = stabilized, sample_info = sample_info)
    expect_false(all(is.na(regressed)))

    de_results <- run(x = regressed, coordinates = coordinates)
    expect_false(all(is.na(de_results)))

    ms_results <- model_search(
        x = regressed, coordinates = coordinates,
        de_results = de_results
    )
    expect_false(all(is.na(ms_results)))

    sp_results <- spatial_patterns(
        x = regressed, coordinates = coordinates,
        de_results = de_results, C = 5L, l = 1.5
    )
    expect_false(all(is.na(sp_results[[1]])))
    expect_false(all(is.na(sp_results[[2]])))

})

test_that("stabilize() warns about NA output", {
    ## Setting a seed to get counts for which stabilize() doesn't work
    set.seed(42)
    counts <- matrix(stats::rpois(ncells * ngenes, lambda = 3),
                     nrow = ngenes, ncol = ncells)

    expect_warning(stabilized <- stabilize(counts))
    expect_true(all(is.na(stabilized)))
})

set.seed(42)
mock <- mockSVG(10, 100, 10)
sample_info <- mock$coordinates
sample_info$total_counts <- colSums(mock$counts)

test_that("stabilize() returns correct output", {
    out <- stabilize(mock$counts)
    expect_equal(nrow(out), nrow(mock$counts))
    expect_true(is.matrix(out))
})

test_that("regress_out() returns correct output", {
    stabilized <- stabilize(mock$counts)
    out <- regress_out(counts = stabilized, sample_info = sample_info)
    expect_equal(nrow(out), nrow(mock$counts))
    expect_equal(nrow(out), nrow(stabilized))
    expect_true(is.matrix(out))
})

test_that("run() returns correct output", {
    stabilized <- stabilize(mock$counts)
    regressed <- regress_out(counts = stabilized, sample_info = sample_info)

    out <- run(x = regressed, coordinates = mock$coordinates)
    expect_equal(nrow(out), nrow(mock$counts))
    expect_true(is.data.frame(out))

    ## Check breaking errors (incompatible dimensions)
    expect_error(run(x = regressed, coordinates = mock$coordinates[1:3, ]))
})

test_that("model_search() and spatial_patterns() return correct output", {
    stabilized <- stabilize(mock$counts)
    regressed <- regress_out(counts = stabilized, sample_info)

    results <- run( x = mock$counts, coordinates = mock$coordinates)
    de_results <- results[results$qval < 0.1, ]

    ## model_search()
    out <- model_search(
        x = regressed, coordinates = mock$coordinates,
        de_results = de_results
    )
    expect_equal(nrow(out), nrow(de_results))
    expect_true(is.data.frame(out))

    C = 1L
    ## spatial_patterns()
    sp <- spatial_patterns(
        x = regressed, coordinates = mock$coordinates,
        de_results = de_results, C = C, l = 1
    )
    pat_res <- sp$pattern_results
    pat <- sp$patterns
    expect_equal(nrow(pat_res), nrow(de_results))
    expect_equal(ncol(pat), C)
    expect_true(all(unique(pat_res$pattern) %in% colnames(pat)))


    ## Check breaking errors (incompatible dimensions)
    expect_error(model_search(
        x = mock$counts, coordinates = mock$coordinates[1:3, ],
        de_results = de_results
    ))
    expect_error(spatial_patterns(
        x = mock$counts, coordinates = mock$coordinates[1:3, ],
        de_results = de_results, C = 2L, l = 1))
})
