#' Wrapper for NaiveDE.stabilize python function
#'
#' @description Wrapped functions from NaiveDE python package
#' [NaiveDE](https://github.com/Teichlab/NaiveDE)
#'
#' @param counts `matrix` or `data.frame` with the expression values.
#' Columns should be samples, and Rows genes
#'
#' @importFrom reticulate import r_to_py
.naiveDE_stabilize <- function(counts) {
    naiveDE <- import("NaiveDE")

    df_py <- r_to_py(as.data.frame(counts))

    stabilized <- naiveDE$stabilize(df_py)

    ## Return stabilized counts as matrix
    as.matrix(stabilized)
}

#' Wrapper for NaiveDE.regress_out python function
#'
#' @param sample_info `data.frame` with samples as rows,
#' 'x', 'y' coordinates and  total raw counts as columns.
#'
#' @param stabilized_counts `matrix` or `data.frame` resulting from
#' naiveDE_stabilize
#'
#' @importFrom reticulate import r_to_py
.naiveDE_regress_out <- function(sample_info, stabilized_counts) {
    naiveDE <- import("NaiveDE")

    sample_info_py <- r_to_py(sample_info)
    df_py <- r_to_py(as.data.frame(stabilized_counts))

    res <- naiveDE$regress_out(sample_info_py, df_py, "np.log(total_counts)")

    as.matrix(res)
}
