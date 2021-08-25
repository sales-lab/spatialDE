## Files originally retrieved from
## https://github.com/Teichlab/SpatialDE/blob/master/Analysis/MouseOB/MOB_sample_info.csv

url <- "https://sales.bio.unipd.it/bulk/b1678ba2e04342053efd5d3bc277d98e65280800ff4da0a1a3fc64e74c3214e8/MOB_sample_info.csv"
MOB_sample_info <- read.csv(url, row.names = 1)

usethis::use_data(MOB_sample_info, overwrite = TRUE)
