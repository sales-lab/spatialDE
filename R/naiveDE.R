#' @description Wrapped functions from NaiveDE python package
#' [NaiveDE](https://github.com/Teichlab/NaiveDE)
#' 
#' @importFrom reticulate import r_to_py
naiveDE_stabilize <- function(df) {
  naiveDE <- import("NaiveDE")
  
  dft <- as.data.frame(t(df))
  dft_py <- r_to_py(dft)
  
  stabilized <- naiveDE$stabilize(dft_py)
  
  tstabilized <- as.data.frame(t(stabilized))
  return(tstabilized)
}

#' @importFrom reticulate import r_to_py
naiveDE_regress_out <- function(sample_info, dfm) {
  naiveDE <- import("NaiveDE")
  
  sample_info_py <- r_to_py(sample_info)
  tdfm_py <- r_to_py(as.data.frame(t(dfm)))
  
  res <- naiveDE$regress_out(sample_info_py, tdfm_py, 'np.log(total_counts)')
  
  tres <- as.data.frame(t(res))
  return(tres)
}
