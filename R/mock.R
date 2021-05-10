#' Generate count matrix for spatially variable genes.
#'
#' @param size An `integer` scalar. Cells will be spatially arranged on a `size
#'   x size` grid. Default: 10, corresponding to 100 cells.
#' @param tot_genes An `integer` scalar. Total number of genes. Default: 1000.
#' @param de_genes An `integer` scaler. The number of spatially variable genes.
#'   Default: 100.
#' @param return_SPE A `logical`, whether to return result as a
#'   \linkS4class{SpatialExperiment}. Default: `FALSE`.
#'
#' @return
#' If `return_SPE = TRUE`, returns a \linkS4class{SpatialExperiment} object.
#'
#' If not, a `list` containing:
#' * `coordinates`: `data.frame` with `x` and `y` columns;
#' * `counts`: `matrix` with generated gene counts.
#'
#' @examples
#' spe <- mockSVG(10, tot_genes = 200, de_genes = 20, return_SPE = TRUE)
#' spe
#'
#' @export
#' @importFrom stats rnbinom runif
mockSVG <- function(size = 10, tot_genes = 1000, de_genes = 100,
                    return_SPE = FALSE) {
    n_cells <- size * size
    coordinates <- data.frame(
        x = rep(seq.int(size), size),
        y = rep(seq.int(size), each = size)
    )

    mu <- 2^runif(tot_genes, -1, 5)
    counts <- matrix(
        rnbinom(tot_genes * n_cells, mu = mu, size = 10),
        nrow = tot_genes
    )
    m <- (size / 2) + 1
    mask <- coordinates$x < m & coordinates$y < m
    counts[seq.int(de_genes), mask] <- counts[seq.int(de_genes), mask] + 20

    ## Set dimnames
    gene_names <- sprintf("Gene_%s", formatC(seq_len(tot_genes), width = 4, flag = 0))
    cell_names <- sprintf("Cell_%s", formatC(seq_len(n_cells), width = 3, flag = 0))
    rownames(coordinates) <- cell_names
    dimnames(counts) <- list(gene_names, cell_names)

    out <- list(coordinates = coordinates, counts = counts)

    if (return_SPE) {
        out <- SpatialExperiment::SpatialExperiment(
            assays = list(counts = counts),
            spatialData = S4Vectors::DataFrame(coordinates),
            spatialCoordsNames = c("x", "y")
        )
    }
    out
}
