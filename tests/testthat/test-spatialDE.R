ncells <- 10
ngenes <- 100

## Note: `stabilize()` returns all NA values for some count inputs
## Setting a seed to get the same counts for which stabilize() works for now
set.seed(1)
counts <- matrix(stats::rpois(ncells * ngenes, lambda = 3),
                 nrow = ngenes, ncol = ncells)
sample_info <- data.frame(total_counts = colSums(counts))
coordinates <- data.frame(x = stats::rnorm(ncells), y = stats::rnorm(ncells))

test_that("Wrapper functions work", {
    stabilized <- stabilize(counts)
    expect_false(all(is.na(stabilized)))

    regressed <- regress_out(sample_info, stabilized)
    expect_false(all(is.na(regressed)))

    de_results <- run(coordinates, regressed)
    expect_false(all(is.na(de_results)))

    ms_results <- model_search(coordinates, regressed, de_results)
    expect_false(all(is.na(ms_results)))

    sp_results <- spatial_patterns(
        coordinates = coordinates, regressed_counts = regressed,
        sres = de_results, C = 5L, l = 1.5
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

test_that("run() returns correct output", {
    out <- run(coordinates = coordinates, regressed_counts = counts)
    expect_equal(nrow(out), nrow(counts))
    expect_true(is.data.frame(out))

    ## Check breaking errors (incompatible dimensions)
    expect_error(run(coordinates = coordinates[1:3, ],
                     regressed_counts = counts))
})

test_that("model_search() and spatial_patterns() return correct output", {
    de_res <- run(coordinates, counts)

    ## model_search()
    out <- model_search(coordinates = coordinates,
                        regressed_counts = counts,
                        de_results = de_res)
    expect_equal(nrow(out), nrow(counts))
    expect_true(is.data.frame(out))

    ## spatial_patterns()
    sp <- spatial_patterns(
        coordinates = coordinates, regressed_counts = counts,
        sres = de_res, C = 1L, l = 1
    )
    pat_res <- sp$pattern_results
    pat <- sp$patterns
    expect_equal(nrow(pat_res), nrow(counts))
    expect_equal(nrow(pat), ncol(counts))


    ## Check breaking errors (incompatible dimensions)
    expect_error(model_search(coordinates = coordinates[1:3, ],
                              regressed_counts = counts,
                              de_results = de_res))
    expect_error(spatial_patterns(coordinates = coordinates[1:3, ],
                              regressed_counts = counts,
                              sres = de_res, C = 2L, l = 1))
})
