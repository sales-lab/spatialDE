#' Plot Fraction Spatial Variance vs Q-value
#'
#' @param results results from SpatialDE.
#' @param ms_results model selection results, should be a data frame with
#'   columns `g` for gene names and `model` for the model selected.
#' @param certain_only only plot results with narrow 95% confidence interval.
#' @param log_x Whether to display x axis in log scale.
#' @param do_label display gene names for statistically significant genes,
#'   default `TRUE`.
#' @param covariate_names names of covariates as a reference, default to `NULL`.
#'
#' @return A `ggplot2` object.
#'
#' @examples
#' set.seed(42)
#' spe <- mockSVG(size = 10, tot_genes = 200, de_genes = 20, return_SPE = TRUE)
#' ## Run spatialDE with S4 integration
#' results <- spatialDE(spe)
#' ## Run model search
#' msearch <- modelSearch(spe, de_results = results, qval_thresh = NULL,
#'   verbose = FALSE)
#'
#' plot <- FSV_sig(results, msearch)
#'
#' @references
#' Svensson, V., Teichmann, S. & Stegle, O. SpatialDE: identification of
#' spatially variable genes. Nat Methods 15, 343–346 (2018).
#' \url{https://doi.org/10.1038/nmeth.4636}
#'
#' [**SpatialDE 1.1.3**](https://pypi.org/project/SpatialDE/1.1.3/): the version
#' of the Python package used under the hood.
#'
#' @author Davide Corso, Milan Malfait, Lambda Moses
#'
#' @import ggplot2
#' @export
FSV_sig <- function(results, ms_results = NULL, certain_only = FALSE,
                    log_x = FALSE, do_label = TRUE, covariate_names = NULL) {
  if (!is.null(ms_results)) {
      results <- merge(results, ms_results[, c("g", "model")],
                       by = "g", all = TRUE, suffixes = c("", "_bic"))
  } else {
      names(results)[names(results) == "model"] <- "model_bic"
  }
  results$FSV95conf <- 2 * sqrt(results$s2_FSV)
  results$conf_categories <- cut(results$FSV95conf, c(0, 0.1, 1, Inf))
  ll <- levels(results$conf_categories)
  results$conf_categories <- factor(results$conf_categories, levels = rev(ll))
  results$is_covariate <- FALSE
  results$color_categories <- ifelse(results$model_bic == "SE", "general",
                                     ifelse(results$model_bic == "PER",
                                            "periodic", "linear"))
  if (!is.null(covariate_names)) {
      results$is_covariate <- results$g %in% covariate_names
  }
  if (certain_only) {
      results <- results[results$conf_categories == "(0,0.1]",]
  }
  FSV <- qval <- color_categories <- conf_categories <- is_covariate <-
      g <- NULL
  colors_use <- scales::hue_pal()(length(unique(results$model_bic)))
  p <- ggplot(results, aes(FSV, qval)) +
      geom_hline(yintercept = 0.05, linetype = 2) +
      scale_color_manual(
          values = colors_use, na.translate = TRUE,
          na.value = "black",
          guide = guide_legend(title = "model")
      ) +
      scale_y_continuous(trans = .reverse_log10()) +
      annotate(geom = "text", x = 0, y = 0.05, label = "0.05")
  if (!is.null(covariate_names)) {
      p <- p +
          scale_shape_manual(
              values = c(16, 4),
              guide = guide_legend(title = "covariate")
          )
  }
  if (!certain_only) {
      p <- p +
          geom_point(aes(
              color = color_categories, size = conf_categories,
              shape = is_covariate
          ), alpha = 0.5) +
          scale_size_manual(
              values = c(0.7, 1.5, 3),
              guide = guide_legend(title = "confidence \n category")
          )
  } else {
      p <- p +
          geom_point(
              aes(color = color_categories, shape = is_covariate),
              alpha = 0.5
          )
  }
  if (log_x) {
      p <- p +
          scale_x_log10()
  }
  if (do_label) {
      gene_label <- results[results$qval < 0.05, c("FSV", "qval", "g")]
      p <- p + ggrepel::geom_label_repel(aes(label = g), data = gene_label)
  }
  p
}

# Inverse log scale
#
# Custom transform for inverted log transformed axis in ggplot2. 
# Return a ggplot2 compatible transformation object
.reverse_log10 <- function() {
    trans <- function(x) -log10(x)
    inv <- function(x) 10^(-x)
    scales::trans_new("reverse_log10", trans, inv,
        scales::log_breaks(base = 10),
        domain = c(1e-100, Inf)
    )
}


#' Plot Spatial Patterns of Multiple Genes
#'
#' @param x A numeric `matrix` of stabilized counts (e.g. resulting from
#' [stabilize()]) where genes are rows and cells are columns.
#'
#'    Alternatively, a \linkS4class{SpatialExperiment} object.
#'
#' @param ... For the generic, arguments to pass to specific methods.
#' @param coordinates A `data.frame` with sample coordinates. Each row is a
#'   sample, the columns with coordinates should be named 'x' and 'y'.
#'
#'   For the *SpatialExperiment* method, coordinates are taken from
#'   `spatialCoords(x)`.
#'
#' @param assay_type A `character` string specifying the assay from `x` to use
#'   as input. Defaults to `"counts"`.
#' @param genes_plot character vector specifying which genes are to be plotted.
#' @param viridis_option This function uses the `viridis` palette to color
#'   cells for gene expression. Four options are available: "magma" (or "A"),
#'   "inferno" (or "B"), "plasma" (or "C"),
#'   "viridis" (or "D", the default option) and "cividis" (or "E").
#' @param ncol Number of columns to arrange the plots.
#' @param point_size Point size of each plot.
#' @param dark_theme Whether dark background should be used; this is helpful to
#'   highlight cells with high expression when using the \code{viridis} palette.
#'
#' @return This function draws a plot for each specified genes
#'
#' @examples
#' ## Mock up a SpatialExperiment object wit 100 cells, 200 genes
#' set.seed(42)
#' spe <- mockSVG(size = 10, tot_genes = 200, de_genes = 10, return_SPE = TRUE)
#'
#' ## Run spatialDE
#' results <- spatialDE(spe)
#'
#' ordered_spe_results <- results[order(results$qval), ]
#' head(ordered_spe_results)
#'
#' plots <- multiGenePlots(spe,
#'     assay_type = "counts",
#'     ordered_spe_results[1:4, ]$g,
#'     point_size = 4,
#'     viridis_option = "D"
#' )
#'
#' @seealso
#' The individual steps performed by this function: [stabilize()],
#' [spatialDE()].
#'
#' For further analysis of the DE results:
#' [model_search()] and [spatial_patterns()].
#'
#' @references
#' Svensson, V., Teichmann, S. & Stegle, O. SpatialDE: identification of
#' spatially variable genes. Nat Methods 15, 343–346 (2018).
#' \url{https://doi.org/10.1038/nmeth.4636}
#'
#' [**SpatialDE 1.1.3**](https://pypi.org/project/SpatialDE/1.1.3/): the version
#' of the Python package used under the hood.
#'
#' @author Davide Corso, Milan Malfait, Lambda Moses
#' @name multiGenePlots
NULL

#' @importFrom gridExtra grid.arrange
#' @importFrom checkmate assert_data_frame assert_names
.multiGenePlots <- function(x, coordinates, genes_plot, viridis_option = "D",
    ncol = 2, point_size = 1, dark_theme = TRUE) {
    assert_data_frame(coordinates, any.missing = FALSE)
    assert_names(colnames(coordinates), identical.to = c("x", "y"))

    stabilized <- x
    y <- NULL
    pls <- lapply(
        seq_along(genes_plot),
        function(i) {
            p <- ggplot(
                data = coordinates,
                aes(
                    x = x, y = y,
                    color = stabilized[genes_plot[i], ]
                )
            ) +
                geom_point(size = point_size) +
                ggtitle(genes_plot[i]) +
                scale_color_viridis_c(option = viridis_option) +
                labs(color = genes_plot[i])

            if (dark_theme) {
                p <- p +
                    theme_dark()
            }
            p
        }
    )
    grid.arrange(grobs = pls, ncol = ncol)
}

#' @import methods
#' @export
#' @rdname multiGenePlots
setGeneric("multiGenePlots", function(x, ...) standardGeneric("multiGenePlots"))

#' @export
#' @rdname multiGenePlots
setMethod("multiGenePlots", "matrix", .multiGenePlots)

#' @export
#' @rdname multiGenePlots
#' @importFrom SummarizedExperiment assay
#' @importFrom SpatialExperiment spatialCoords spatialCoordsNames
#' @importFrom basilisk basiliskStart basiliskRun
setMethod("multiGenePlots", "SpatialExperiment",
    function(x, assay_type = "counts", genes_plot, viridis_option = "D",
             ncol = 2, point_size = 1, dark_theme = TRUE) {
        ## Rename spatialCoords columns to "x", "y"
        spatialCoordsNames(x) <- c("x", "y")
        coordinates <- as.data.frame(spatialCoords(x))

        proc <- basiliskStart(spatialDE_env, testload="scipy.optimize")
        # Stabilize
        counts <- assay(x, assay_type)
        
        assert_matrix(counts, any.missing = FALSE)
        .naiveDE_stabilize(proc, counts)
        stabilized <- basiliskRun(proc, function(store) {
          as.matrix(store$stabilized)
        }, persist=TRUE)

        .multiGenePlots(
            x = stabilized,
            coordinates = coordinates,
            genes_plot = genes_plot,
            viridis_option = viridis_option,
            ncol = ncol,
            point_size = point_size,
            dark_theme = dark_theme
        )
    }
)
