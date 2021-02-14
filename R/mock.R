#' Generate count matrix for spatially variable genes.
#'
#' @param size `int` genes will be spatially arranged on a size x size grid
#' @param tot_genes `int` total number of genes
#' @param de_genes `int` number of spatially variable genes
#'
#' @return `list` containing:
#' * `coordinates`: `data.frame` with `x` and `y` columns;
#' * `counts`: `matrix` with generated gene counts.
#' @export
#' @importFrom stats rnbinom runif
mockSVG <- function(size, tot_genes, de_genes) {
    coordinates <- data.frame(
        x = rep(seq.int(size), size),
        y = rep(seq.int(size), each = size)
    )

    mu <- 2^runif(tot_genes, -1, 5)
    counts <- matrix(
        rnbinom(tot_genes * size * size, mu = mu, size = 10),
        nrow = tot_genes
    )
    m <- (size / 2) + 1
    mask <- coordinates$x < m & coordinates$y < m
    counts[seq.int(de_genes), mask] <- counts[seq.int(de_genes), mask] + 20

    list(coordinates = coordinates, counts = counts)
}
