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

    de_results <- run(x = regressed, coordinates = coordinates)
    expect_false(all(is.na(de_results)))
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

set.seed(42)
mock <- mockSVG(10, 100, 10)
sample_info <- mock$coordinates
sample_info$total_counts <- colSums(mock$counts)

stabilized <- stabilize(mock$counts)
regressed <- regress_out(counts = stabilized, sample_info = sample_info)
de_results <- run(x = regressed, coordinates = mock$coordinates)

test_that("stabilize() returns correct output", {
    expect_equal(nrow(stabilized), nrow(mock$counts))
    expect_true(is.matrix(stabilized))
})

test_that("regress_out() returns correct output", {
    expect_equal(nrow(regressed), nrow(mock$counts))
    expect_equal(nrow(regressed), nrow(stabilized))
    expect_true(is.matrix(regressed))
})

test_that("run() returns correct output", {
    expect_equal(nrow(de_results), nrow(mock$counts))
    expect_true(is.data.frame(de_results))

    ## Check breaking errors (incompatible dimensions)
    expect_error(run(x = regressed, coordinates = mock$coordinates[1:3, ]))
})
