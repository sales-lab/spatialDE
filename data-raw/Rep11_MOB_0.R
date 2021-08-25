library(Matrix)

## File originally retrieved from
## https://github.com/Teichlab/SpatialDE/blob/master/Analysis/MouseOB/data/Rep11_MOB_0.csv

url <- "https://sales.bio.unipd.it/bulk/b1678ba2e04342053efd5d3bc277d98e65280800ff4da0a1a3fc64e74c3214e8/Rep11_MOB_0.csv"
Rep11_MOB_0 <- read.csv(url, header = TRUE, row.names = 1)

## Convert to sparse matrix and transpose to get genes as rows
Rep11_MOB_0 <- as.matrix(Rep11_MOB_0)
Rep11_MOB_0 <- t(Rep11_MOB_0)

usethis::use_data(Rep11_MOB_0, overwrite = TRUE)
