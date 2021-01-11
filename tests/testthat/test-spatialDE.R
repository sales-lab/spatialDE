ncells <- 10
ngenes <- 100

## FIXME: `stabilize()` returns all NA values for some count inputs
## Setting a seed to get the same counts for which stabilize() works for now
set.seed(1)
counts <- matrix(rpois(ncells * ngenes, lambda = 3),
                 nrow = ngenes, ncol = ncells)
sample_info <- data.frame(total_counts = colSums(counts))
coordinates <- data.frame(x = rnorm(ncells), y = rnorm(ncells))

test_that("Wrapper functions work", {
    stabilized <- stabilize(counts)
    expect_false(all(is.na(stabilized)))

    regressed <- regress_out(sample_info, stabilized)
    expect_false(all(is.na(regressed)))

    de_results <- run_spatialDE(coordinates, regressed)
    expect_false(all(is.na(de_results)))

    ms_results <- run_model_search(coordinates, regressed, de_results)
    expect_false(all(is.na(ms_results)))
})
