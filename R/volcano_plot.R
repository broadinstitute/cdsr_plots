#' Make volcano plot
#'
#' @param data data frame containing stats
#' @param effect_var variable name for effect size (x-axis)
#' @param p_var p value - variable name for y-axis
#' @param q_var q value--when specified, defaults to highlighting points that pass q_thresh
#' @param q_thresh q value threshold for highlighting points
#' @param label_var column containing labels for data points
#' @param n_labeled if rank_by = 'effect', n_labeled points per side; if rank_by = 'pval', n_labeled points in total
#' @param label_size size of points
#' @param label_bool logical column indicating which points should be labeled
#' @param rank_by data type used to rank data points when labeling
#' @param ggrepel_type choose 'text' or 'label' to determine whether ggrepel's geom_text_repel or geom_label_repel should be used
#' @param color_var logical/categorical column for coloring points (if a factor, points will be layered according to the levels)
#' @param color_values vector assigning categories from color_var to custom colors
#'
#' @importFrom magrittr "%>%"
#' @importFrom magrittr "%<>%"
#'
#' @examples
#' ceres <- taigr::load.from.taiga(data.name='avana-internal-19q1-7643', data.version=1, data.file='gene_effect')
#' sample_info <- taigr::load.from.taiga(data.name='avana-internal-19q1-7643', data.version=1, data.file='sample_info')
#'
#' groups <- (sample_info$primary_tissue == 'skin') %>% magrittr::set_names(sample_info$DepMap_ID)
#' limma_res <- cdsr::run_lm_stats_limma(ceres[names(groups),], groups) %>% dplyr::arrange(desc(abs(EffectSize)))
#'
#' limma_res %>%
#'     dplyr::mutate(fdr = p.adjust(p.value, method = 'fdr')) %>%
#'     cdsr::make_volcano('EffectSize', 'p.value', q_var = 'fdr')
#' @export
#'
make_volcano <- function(data, effect_var, p_var, q_var = NULL, q_thresh = 0.1, label_var = NULL,
                         n_labeled = 10, label_size = 3, label_bool = NULL, rank_by = c('effect', 'pval'),
                         ggrepel_type = c('text', 'label'), color_var = NULL, color_values = NULL) {
  library(ggplot2)
  # set label for colors in the legend
  guide_title <- color_var
  # log 10 transform the p values
  transformed_p_var <- paste0('-log10(', p_var, ')')
  data[[transformed_p_var]] <- -log10(data[[p_var]])

  # if user has specified q values but no variable to color by, color points that pass q_thresh
  if (is.null(color_var) & !is.null(q_var)) {
    color_var <- 'internal_sig'
    data[[color_var]] <- data[[q_var]] < q_thresh
    guide_title <- sprintf('FDR < %.3f', q_thresh)
  }

  if (is.null(color_var)) {
    p <- data %>%
      ggplot(aes_string(effect_var, transformed_p_var)) +
      geom_point(color = '#333333', alpha = 0.7)
  } else {
    if (is.null(color_values)) { # user has not specified exact colors, create default color schemes
      if (is.logical(data[[color_var]])) color_values <- c(`TRUE` = '#BF2026', `FALSE` = '#4d4d4d')
      else color_values <- RColorBrewer::brewer.pal(length(unique(data[[color_var]])), 'Dark2')
    }
    # plot layers one by one
    if (is.factor(data[[color_var]])) {
      layering <- levels(data[[color_var]])
    } else {
      layering <- sort(unique(as.character(data[[color_var]])))
    }
    p <- ggplot(data, aes_string(effect_var, transformed_p_var, color = color_var))
    for (cur_layer in layering) p <- p + geom_point(data = data[data[[color_var]] == cur_layer,], alpha = 0.7)
    p <- p +
      scale_color_manual(values = color_values) +
      ggplot2::guides(color = guide_legend(title = guide_title))
  }

  if (!is.null(label_var)) { # user has specified column to label points
    if (is.null(label_bool)) { # default to labeling top 10 data points on each side by effect size
      label_bool <- 'to_label' # define points to label with this column
      if (rank_by[1] == 'effect') {
        left <- rank(data[[effect_var]], ties.method = 'random')
        right <- rank(-data[[effect_var]], ties.method = 'random')
        data[[label_bool]] <- (left <= n_labeled) | (right <= n_labeled)
      } else if (rank_by[1] == 'pval') {
        data[[label_bool]] <- rank(data[[p_var]], ties.method = 'random') <= n_labeled
      }
    }
    if (ggrepel_type[1] == 'text') p <- p +
        ggrepel::geom_text_repel(data = data[data[[label_bool]],], aes_string(label = label_var), size = label_size, show.legend = F)
    else if (ggrepel_type[1] == 'label') p <- p +
        ggrepel::geom_label_repel(data = data[data[[label_bool]],], aes_string(label = label_var), size = label_size, show.legend = F)
  }
  return(p)
}


