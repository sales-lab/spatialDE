ncells <- 10
ngenes <- 100

## Note: `stabilize()` returns all NA values for some count inputs
## Setting a seed to get the same counts for which stabilize() works for now
set.seed(1)
counts <- matrix(stats::rpois(ncells * ngenes, lambda = 3),
    nrow = ngenes, ncol = ncells
)
sample_info <- data.frame(
    x = stats::rnorm(ncells), y = stats::rnorm(ncells),
    total_counts = colSums(counts)
)
coordinates <- sample_info[, c("x", "y")]

test_that("Wrapper functions work", {
    stabilized <- stabilize(counts = counts)
    expect_false(all(is.na(stabilized)))

    regressed <- regress_out(counts = stabilized, sample_info = sample_info)
    expect_false(all(is.na(regressed)))
})

test_that("stabilize() warns about NA output", {
    ## Setting a seed to get counts for which stabilize() doesn't work
    set.seed(42)
    counts <- matrix(stats::rpois(ncells * ngenes, lambda = 3),
        nrow = ngenes, ncol = ncells
    )

    expect_warning(stabilized <- stabilize(counts))
    expect_true(all(is.na(stabilized)))
})
