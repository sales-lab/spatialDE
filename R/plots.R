
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
#' @import ggplot2 forcats
#' @importFrom dplyr mutate filter select full_join case_when
#' @importFrom magrittr %>%
#' @export
FSV_sig <- function(results, ms_results = NULL, certain_only = FALSE, 
                    log_x = FALSE, do_label = TRUE, covariate_names = NULL) {
  if (!is.null(ms_results)) {
    results <- results %>%
      full_join(ms_results[,c("g", "model")], by = "g", suffix = c("", "_bic"))
  } else {
    results <- results %>%
      rename(model_bic = model)
  }
  results <- results %>%
    mutate(FSV95conf = 2 * sqrt(s2_FSV),
           conf_categories = fct_rev(cut(FSV95conf, c(0, 0.1, 1, Inf))),
           is_covariate = FALSE,
           # More user friendly model labels
           color_categories = case_when(model_bic == "SE" ~ "general",
                                        model_bic == "PER" ~ "periodic",
                                        model_bic == "linear" ~ "linear"))
  if (!is.null(covariate_names)) {
    results <- results %>%
      mutate(is_covariate = g %in% covariate_names)
  }
  if (certain_only) {
    results <- results %>%
      filter(conf_categories == "(0,0.1]")
  }
  colors_use <- scales::hue_pal()(length(unique(results$model_bic)))
  p <- ggplot(results, aes(FSV, qval)) +
    geom_hline(yintercept = 0.05, linetype = 2) +
    scale_color_manual(values = colors_use, na.translate = TRUE, 
                       na.value = "black", 
                       guide = guide_legend(title = "model")) +
    scale_y_continuous(trans = .reverse_log10()) +
    annotate(geom = "text", x = 0, y = 0.05, label = "0.05")
  if (!is.null(covariate_names)) {
    p <- p +
      scale_shape_manual(values = c(16,4), 
                         guide = guide_legend(title = "covariate"))
  }
  if (!certain_only) {
    p <- p +
      geom_point(aes(color = color_categories, size = conf_categories, 
                     shape = is_covariate), alpha = 0.5) +
      scale_size_manual(values = c(0.7, 1.5, 3),
                        guide = guide_legend(title = "confidence \n category"))
  } else {
    p <- p +
      geom_point(aes(color = color_categories, shape = is_covariate),
                 alpha = 0.5)
  }
  if (log_x) {
    p <- p +
      scale_x_log10()
  }
  if (do_label) {
    gene_label <- results %>%
      filter(qval < 0.05) %>%
      select(FSV, qval, g)
    p <- p + ggrepel::geom_label_repel(aes(label = g), data = gene_label)
  }
  p
}

#' Inverse log scale
#'
#' Custom transform for inverted log transformed axis in ggplot2.
#' @importFrom scales trans_new log_breaks
.reverse_log10 <- function() {
  trans <- function(x) -log10(x)
  inv <- function(x) 10^(-x)
  scales::trans_new("reverse_log10", trans, inv, scales::log_breaks(base = 10),
                    domain = c(1e-100, Inf))
}
