#' Stabilize variance of counts
#'
#' Stabilize variance of negative binomial data using Anscombe's approximation.
#' This function is a wrapper for `stabilize` from the
#' [NaiveDE](https://github.com/Teichlab/NaiveDE) Python package.
#'
#' @param counts `matrix` with expression values for samples in
#' columns and genes in rows.
#'
#' @return `matrix` of variance stabilized counts.
#'
#' @examples
#' set.seed(42)
#' mock <- mockSVG(10, 1000, 10)
#' stabilized <- stabilize(mock$counts)
#' @export
#' @importFrom checkmate assert_matrix test_matrix
#' @importFrom basilisk basiliskStart basiliskRun
stabilize <- function(counts) {
    assert_matrix(counts, any.missing = FALSE)
    proc <- basiliskStart(spatialDE_env, testload="scipy.optimize")
    
    .naiveDE_stabilize(proc, counts)
    
    out <- basiliskRun(proc, function(store) {
        as.matrix(store$stabilized)
    }, persist=TRUE)

    if (!test_matrix(out, all.missing = FALSE)) {
        warning("Warning: Stabilized values are all NA.\n")
    }

    out
}


#' @importFrom reticulate import r_to_py
#' @importFrom basilisk basiliskRun
.naiveDE_stabilize <- function(proc, counts) {
    basiliskRun(proc, function(counts, store) {
        naiveDE <- import("NaiveDE")
        
        ## NaiveDE.stabilize requires data.frame input to work
        df_py <- r_to_py(as.data.frame(counts))
        stabilized <- naiveDE$stabilize(df_py)
        
        store$stabilized <- as.matrix(stabilized)
        invisible(NULL)
    }, counts=counts, persist=TRUE)
}


#' Regress out library size effect
#'
#' Regresses out the effect of library size.
#' This function is a wrapper for `regress_out` from the
#' [NaiveDE](https://github.com/Teichlab/NaiveDE) Python package.
#'
#' @param counts `matrix` of variance stabilized
#' counts, e.g. resulting from [stabilize()].
#' @param sample_info `data.frame` with samples as rows and at least a column
#' with `total_counts`.
#'
#' @return `matrix` of normalized counts.
#'
#' @examples
#' set.seed(42)
#' mock <- mockSVG(10, 1000, 10)
#' stabilized <- stabilize(mock$counts)
#' sample_info <- mock$coordinates
#' sample_info$total_counts <- colSums(mock$counts)
#' regressed <- regress_out(counts = stabilized, sample_info = sample_info)
#' @export
#' @importFrom checkmate assert_data_frame assert_names assert_matrix
#' @importFrom basilisk basiliskStart basiliskRun
regress_out <- function(counts, sample_info) {
    assert_data_frame(sample_info, any.missing = FALSE)
    assert_names(colnames(sample_info),
        must.include = "total_counts"
    )
    assert_matrix(counts, any.missing = FALSE)
    
    proc <- basiliskStart(spatialDE_env, testload="scipy.optimize")
    .naiveDE_regress_out(proc, counts, sample_info)
    
    out <- basiliskRun(proc, function(store) {
        as.matrix(store$regressed)
    }, persist=TRUE)
    
    out
}


#' @importFrom reticulate import r_to_py
#' @importFrom basilisk basiliskRun
.naiveDE_regress_out <- function(proc, counts, sample_info) {
    basiliskRun(proc, function(counts, sample_info, store) {
        naiveDE <- import("NaiveDE")
        
        sample_info_py <- r_to_py(sample_info)
        df_py <- r_to_py(as.data.frame(counts))
        
        regressed <- naiveDE$regress_out(sample_info_py, 
                                         df_py, 
                                         "np.log(total_counts)")
        
        store$regressed <- as.matrix(regressed)
        
        invisible(NULL)
    }, counts=counts, sample_info=sample_info, persist=TRUE)
}
