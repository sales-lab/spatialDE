#' Stabilize variance of counts
#'
#' Stabilize variance of negative binomial data using Anscombe's approximation.
#' This function is a wrapper for `stabilize` from the
#' [NaiveDE](https://github.com/Teichlab/NaiveDE) Python package.
#'
#' @param counts `matrix` or `data.frame` with expression values for samples in
#' columns and genes in rows.
#'
#' @return `matrix` of variance stabilized counts.
#'
#' @examples
#' ncells <- 100
#' ngenes <- 1000
#' counts <- matrix(rpois(ncells * ngenes, lambda = 3),
#'                  nrow = ngenes, ncol = ncells)
#' stabilized <- stabilize(counts)
#'
#' @export
stabilize <- function(counts) {
    out <- basilisk::basiliskRun(
        env = spatialDE_env,
        fun = .naiveDE_stabilize,
        counts = counts
    )

    if (all(is.na(out))) {
        warning("Stabilized values are all NA.")
    }

    out
}


#' @importFrom reticulate import r_to_py
.naiveDE_stabilize <- function(counts) {
    naiveDE <- import("NaiveDE")

    ## NaiveDE.stabilize requires data.frame input to work
    df_py <- r_to_py(as.data.frame(counts))

    stabilized <- naiveDE$stabilize(df_py)

    ## Return stabilized counts as matrix
    as.matrix(stabilized)
}



#' Regress out library size effect
#'
#' Regresses out the effect of library size.
#' This function is a wrapper for `regress_out` from the
#' [NaiveDE](https://github.com/Teichlab/NaiveDE) Python package.
#'
#' @param sample_info `data.frame` with samples as rows and at least a column
#' with `total_counts`.
#'
#' @param stabilized_counts `matrix` or `data.frame` of variance stabilized
#' counts, e.g. resulting from [stabilize()].
#'
#' @return `matrix` of normalized counts.
#'
#' @examples
#' ncells <- 100
#' ngenes <- 1000
#' counts <- matrix(rpois(ncells * ngenes, lambda = 3),
#'                  nrow = ngenes, ncol = ncells)
#'
#' ## Provide total counts for library size normalization
#' sample_info <- data.frame(total_counts = colSums(counts))
#'
#' stabilized <- stabilize(counts)
#' regressed <- regress_out(sample_info, stabilized)
#'
#' @export
regress_out <- function(sample_info, stabilized_counts) {
    out <- basilisk::basiliskRun(
        env = spatialDE_env,
        fun = .naiveDE_regress_out,
        sample_info = sample_info, stabilized_counts = stabilized_counts
    )
    out
}


#' @importFrom reticulate import r_to_py
.naiveDE_regress_out <- function(sample_info, stabilized_counts) {
    naiveDE <- import("NaiveDE")

    # TODO: add check whether `sample-info` contains `total_counts` column

    sample_info_py <- r_to_py(sample_info)

    ## NaiveDE.regress_out requires data.frame input to work
    df_py <- r_to_py(as.data.frame(stabilized_counts))


    # FIXME: the call below only uses the `total_counts` column from sample_info
    # to fit the linear model. So it will probably be safer if we just create
    # sample_info using the provided counts.

    res <- naiveDE$regress_out(sample_info_py, df_py, "np.log(total_counts)")

    as.matrix(res)
}
