#' Make dot plot with gene set enrichment results
#'
#' @param data a table containing gene set enrichment results.
#' @param enrich_var name of the column containing enrichment values. Defaults to "NES" or "odds_ratio".
#' @param size_var name of the column containing sizes. Defaults to "size" or "overlap_size".
#' @param p_var name of the column containing p-values. Defaults to "p_value".
#' @param dir direction(s) to plot. Options are ("both","pos","neg").
#' @param color_by variable to color points by. Options are ("pval","dir","enrich").
#' @param x_by variable to plot on x-axis. Options are ("pval","enrich").
#' @param n_shown number of terms to show.
#' @param sig_only whether to only show significantly enriched terms.
#' @param term_var name of the column containing the term. Defaults to "term".
#' @param remove_collection_name whether to remover the gene set collection name e.g. HALLMARK.
#' @param term_words_per_line number of words per line in term.
#' @param term_max_lines maximum number of lines the term can use.
#' @param theme ggplot theme applied to the plot.
#'
#' @importFrom magrittr "%>%"
#' @importFrom magrittr "%<>%"
#'
#' @export make_gsea_dot
make_gsea_dot <- function(data, enrich_var = NULL, size_var = NULL, p_var = "p_value",dir = "both",
                          color_by = "pval", x_by = "enrich", n_shown = 10, sig_only = FALSE, term_var = "term",
                          remove_collection_name = FALSE, term_words_per_line = 2, term_max_lines = 4,
                          theme = NULL) {
  library(tidyverse)
  library(cowplot)
  #Set defaults
  if(is.null(enrich_var)) {
    if("NES" %in% names(data)) {enrich_var <- "NES"}
    else if ("odds_ratio" %in% names(data)) {enrich_var <- "odds_ratio"}
    else{stop("enrich_var must specified since neither NES nor odds_ratio is a column in data")}
  }
  if(is.null(size_var)) {
    if("overlap_size" %in% names(data)) {size_var <- "overlap_size"}
    else if ("size" %in% names(data)) {size_var <- "size"}
    else{stop("size_var must specified since neither size nor overlap_size is a column in data")}
  }
  #Check inputs
  if(!enrich_var %in% names(data)){stop("enrich_var must be a column in data")}
  if(!size_var %in% names(data)){stop("size_var must be a column in data")}
  if(!p_var %in% names(data)){stop("p_var must be a column in data")}
  if(!dir %in% c("pos","neg","both")) {stop("dir must be in c('pos','neg','both')")}
  if(!color_by %in% c("pval","dir","enrich")) {stop("color_by must be in c('pval','dir','enrich')")}
  if(!x_by %in% c("pval","enrich")) {stop("x_by must be either 'pval' or 'enrich'")}
  if(!term_var %in% names(data)){stop("term_var must be a column in data")}
  if(!"direction" %in% names(data)) {
    if ("NES" %in% names(data)) {
      data[["direction"]] <- ifelse(data$NES > 0,"pos","neg")
    } else if (dir != "both") {
      data[["direction"]] <- dir
    }
    else {
      stop("direction must be a column in data if dir is 'both'")
    }
  }
  #Process inputs
  transformed_p_var <- str_c("-log10(",p_var,")")
  data[transformed_p_var] <- -log10(data[[p_var]])
  data <- data %>% arrange(-!!as.symbol(transformed_p_var))
  if(color_by == "pval") {color_var <- transformed_p_var}
  else if(color_by == "dir") {color_var <- "direction"}
  else {
    if (enrich_var %in% c("NES","ES")) {
      color_var <- str_c("abs(",enrich_var,")")
      data[[color_var]] <- abs(data[[enrich_var]])
    }
    else (color_var <- enrich_var)
  }
  x_var <- ifelse(x_by == "enrich",enrich_var,transformed_p_var)
  #Process terms
  term <- data[[term_var]]
  term %<>% str_replace_all("_"," ")
  if(remove_collection_name) {term %<>% str_remove_all("REACTOME |BIOCARTA |HALLMARK |KEGG |PID ")}
  term %<>% str_replace_all(str_c("(",strrep("[^ ]+ ",term_words_per_line),")"), "\\1\n")
  term_split <- term %>% str_split_fixed(" \n",term_max_lines + 1)
  term <- term_split[,1]
  for (i in 2:term_max_lines) {
    term %<>% str_c(ifelse(term_split[,i] == "","","\n"),term_split[,i])
  }
  term %<>% str_c(ifelse(term_split[,term_max_lines + 1] == "","","..."))
  data["term"] <- term
  #Positive or Negitive
  if (dir != "both") {
    data <- data %>% filter(direction == dir) %>% head(n_shown)
    if (sig_only) {data <- filter(data,!!as.symbol(p_var) < .05)}
    if (nrow(data) == 0){
      return(ggdraw() + draw_text("No terms are significantly enriched", .5,.6))
    }
    max_x <- max(abs(data[[x_var]]))
    p <- data %>% mutate(term = factor(term,levels = rev(term))) %>%
      ggplot(aes_string(x = x_var,y = "term",size = size_var,color = color_var)) +
      geom_point() +  ylab("") + lims(x = c(0,ifelse(data[[1,x_var]] < 0,-max_x,max_x)))
    if (color_var == "direction") {
      p <- p + scale_colour_manual(values = c("pos" = "red", "neg" = "blue")) + guides(color = F)
    } else {
      p <- p + scale_color_gradient(low = "blue", high = "red",limits = c(0,max(data[[color_var]])))
    }
    p <- p + theme
    #Both
  } else {
    data_neg <- data %>% filter(direction == "neg") %>% head(n_shown)
    data_pos <- data %>% filter(direction == "pos") %>% head(n_shown)
    if (sig_only) {
      data_neg <- filter(data_neg,!!as.symbol(p_var) < .05)
      data_pos <- filter(data_pos,!!as.symbol(p_var) < .05)
    }
    data <- rbind(data_neg,data_pos)
    if (nrow(data) == 0){
      return(ggdraw() + draw_text("No terms are significantly enriched", .5,.6))
    }
    max_x <- max(abs(data[[x_var]]))
    max_size <- max(data[[size_var]])
    max_color <- max(data[[color_var]])
    if (nrow(data_neg) == 0) {
      p_neg <- draw_text("No terms are significantly\n enriched in the negative\n direction.", .25,.6)
    } else {
      p_neg <- data_neg %>% mutate(term = factor(term,levels = rev(term))) %>%
        ggplot(aes_string(x_var,"term",size = size_var,color = color_var)) +
        geom_point() +  ylab("") +
        lims(size = c(0,max_size),x = c(ifelse(data_neg[[1,x_var]] < 0,-max_x,max_x),0))
      if (color_var == "direction") {
        p_neg <- p_neg + scale_colour_manual(values = c("pos" = "red", "neg" = "blue")) + guides(color = F)
      } else {
        p_neg <- p_neg + scale_color_gradient(low = "blue", high = "red",limits = c(0,max_color))
      }
      p_neg <- p_neg + ggtitle("Negative") + theme(text = element_text(size=9),plot.margin = unit(c(5,2,1,1), "pt"),
                                                   legend.position="bottom",plot.title = element_text(hjust = 0.5)) + theme
      legend <- cowplot::get_legend(p_neg)
      p_neg <- draw_plot(p_neg + theme(legend.position="none"), 0, .15, .5, .85)
    }
    if (nrow(data_pos) == 0) {
      p_pos <- draw_text("No terms are significantly\n enriched in the positive\n direction.", .75, .6)
    } else {
      p_pos <- data_pos %>% mutate(term = factor(term,levels = rev(term))) %>%
        ggplot(aes_string(x_var,"term",size = size_var,color = color_var)) +
        geom_point() +  ylab("") + scale_y_discrete(position = "right") +
        lims(size = c(0,max_size), x = c(0,max_x))
      if (color_var == "direction") {
        p_pos <- p_pos + scale_colour_manual(values = c("pos" = "red", "neg" = "blue")) + guides(color = F)
      } else {
        p_pos <- p_pos + scale_color_gradient(low = "blue", high = "red",limits = c(0,max_color))
      }
      p_pos <- p_pos + ggtitle("Positive") +theme(text = element_text(size=9), plot.margin = unit(c(5,1,1,2), "pt"),
                                                  legend.position="bottom",plot.title = element_text(hjust = 0.5)) + theme
      legend <- cowplot::get_legend(p_pos)
      p_pos <- draw_plot(p_pos + theme(legend.position="none"), .5, .15, .5, .85)
    }
    p <-ggdraw() + p_neg + p_pos + draw_grob(legend, 0, 0, 1, .15)
  }
  return(p)
}

#' Make bar plot with gene set enrichment results
#'
#' @param data a table containing gene set enrichment results.
#' @param enrich_var name of the column containing enrichment values. Defaults to "NES" or "odds_ratio".
#' @param p_var name of the column containing p-values. Defaults to "p_value".
#' @param dir direction(s) to plot. Options are ("both","pos","neg").
#' @param color_by variable to color points by. Options are ("pval","dir","enrich").
#' @param x_by variable to plot on x-axis. Options are ("pval","enrich").
#' @param n_shown number of terms to show.
#' @param sig_only whether to only show significantly enriched terms.
#' @param term_var name of the column containing the term. Defaults to "term".
#' @param remove_collection_name whether to remover the gene set collection name e.g. HALLMARK.
#' @param term_words_per_line number of words per line in term.
#' @param term_max_lines maximum number of lines the term can use.
#' @param theme ggplot theme applied to the plot.
#'
#' @importFrom magrittr "%>%"
#' @importFrom magrittr "%<>%"
#'
#' @export make_gsea_bar
make_gsea_bar <- function(data, enrich_var = NULL, p_var = "p_value",dir = "both", color_by = "pval",
                          x_by = "enrich", n_shown = 10, sig_only = FALSE, term_var = "term",
                          remove_collection_name = FALSE, term_words_per_line = 2, term_max_lines = 4,
                          theme = NULL) {
  library(tidyverse)
  library(cowplot)
  #Set defaults
  if(is.null(enrich_var)) {
    if("NES" %in% names(data)) {enrich_var <- "NES"}
    else if ("odds_ratio" %in% names(data)) {enrich_var <- "odds_ratio"}
    else{stop("enrich_var must specified since neither NES nor odds_ratio is a column in data")}
  }
  #Check inputs
  if(!enrich_var %in% names(data)){stop("enrich_var must be a column in data")}
  if(!p_var %in% names(data)){stop("p_var must be a column in data")}
  if(!dir %in% c("pos","neg","both")) {stop("dir must be in c('pos','neg','both')")}
  if(!color_by %in% c("pval","dir","enrich")) {stop("color_by must be in c('pval','dir','enrich')")}
  if(!x_by %in% c("pval","enrich")) {stop("x_by must be either 'pval' or 'enrich'")}
  if(!term_var %in% names(data)){stop("term_var must be a column in data")}
  if(!"direction" %in% names(data)) {
    if ("NES" %in% names(data)) {
      data[["direction"]] <- ifelse(data$NES > 0,"pos","neg")
    } else if (dir != "both") {
      data[["direction"]] <- dir
    }
    else {
      stop("direction must be a column in data if dir is 'both'")
    }
  }
  #Process inputs
  transformed_p_var <- str_c("-log10(",p_var,")")
  data[transformed_p_var] <- -log10(data[[p_var]])
  data <- data %>% arrange(-!!as.symbol(transformed_p_var))
  if(color_by == "pval") {color_var <- transformed_p_var}
  else if(color_by == "dir") {color_var <- "direction"}
  else {
    if (enrich_var %in% c("NES","ES")) {
      color_var <- str_c("abs(",enrich_var,")")
      data[[color_var]] <- abs(data[[enrich_var]])
    }
    else (color_var <- enrich_var)
  }
  x_var <- ifelse(x_by == "enrich",enrich_var,transformed_p_var)
  #Process terms
  term <- data[[term_var]]
  term %<>% str_replace_all("_"," ")
  if(remove_collection_name) {term %<>% str_remove_all("REACTOME |BIOCARTA |HALLMARK |KEGG |PID ")}
  term %<>% str_replace_all(str_c("(",strrep("[^ ]+ ",term_words_per_line),")"), "\\1\n")
  term_split <- term %>% str_split_fixed(" \n",term_max_lines + 1)
  term <- term_split[,1]
  for (i in 2:term_max_lines) {
    term %<>% str_c(ifelse(term_split[,i] == "","","\n"),term_split[,i])
  }
  term %<>% str_c(ifelse(term_split[,term_max_lines + 1] == "","","..."))
  data["term"] <- term
  #Positive or Negitive
  if (dir != "both") {
    data <- data %>% filter(direction == dir) %>% head(n_shown)
    if (sig_only) {data <- filter(data,!!as.symbol(p_var) < .05)}
    if (nrow(data) == 0){
      return(ggdraw() + draw_text("No terms are significantly enriched", .5,.6))
    }
    max_x <- max(abs(data[[x_var]]))
    p <- data %>% mutate(term = factor(term,levels = rev(term))) %>%
      ggplot(aes_string(x = x_var,y = "term",fill = color_var)) +
      geom_col() +  ylab("") + lims(x = c(0,ifelse(data[[1,x_var]] < 0,-max_x,max_x)))
    if (color_var == "direction") {
      p <- p + scale_colour_manual(values = c("pos" = "red", "neg" = "blue")) + guides(fill = F)
    } else {
      p <- p + scale_color_gradient(low = "blue", high = "red",limits = c(0,max(data[[color_var]])))
    }
    p <- p + theme
    #Both
  } else {
    data_neg <- data %>% filter(direction == "neg") %>% head(n_shown)
    data_pos <- data %>% filter(direction == "pos") %>% head(n_shown)
    if (sig_only) {
      data_neg <- filter(data_neg,!!as.symbol(p_var) < .05)
      data_pos <- filter(data_pos,!!as.symbol(p_var) < .05)
    }
    data <- rbind(data_neg,data_pos)
    if (nrow(data) == 0){
      return(ggdraw() + draw_text("No terms are significantly enriched", .5,.6))
    }
    max_x <- max(abs(data[[x_var]]))
    max_color <- max(data[[color_var]])
    if (nrow(data_neg) == 0) {
      p_neg <- draw_text("No terms are significantly\n enriched in the negative\n direction.", .25,.6)
    } else {
      p_neg <- data_neg %>% mutate(term = factor(term,levels = rev(term))) %>%
        ggplot(aes_string(x_var,"term",fill = color_var)) +
        geom_col() +  ylab("") +
        lims(x = c(ifelse(data_neg[[1,x_var]] < 0,-max_x,max_x),0))
      if (color_var == "direction") {
        p_neg <- p_neg + scale_fill_manual(values = c("pos" = "red", "neg" = "blue")) + guides(fill = F)
      } else {
        p_neg <- p_neg + scale_fill_gradient(low = "blue", high = "red",limits = c(0,max_color))
      }
      p_neg <- p_neg + ggtitle("Negative") + theme(text = element_text(size=9),plot.margin = unit(c(5,2,1,1), "pt"),
                                                   legend.position="bottom",plot.title = element_text(hjust = 0.5)) + theme
      legend <- cowplot::get_legend(p_neg)
      p_neg <- draw_plot(p_neg + theme(legend.position="none"), 0, .15, .5, .85)
    }
    if (nrow(data_pos) == 0) {
      p_pos <- draw_text("No terms are significantly\n enriched in the positive\n direction.", .75, .6)
    } else {
      p_pos <- data_pos %>% mutate(term = factor(term,levels = rev(term))) %>%
        ggplot(aes_string(x_var,"term",fill = color_var)) +
        geom_col() +  ylab("") + scale_y_discrete(position = "right") +
        lims(x = c(0,max_x))
      if (color_var == "direction") {
        p_pos <- p_pos + scale_fill_manual(values = c("pos" = "red", "neg" = "blue")) + guides(fill = F)
      } else {
        p_pos <- p_pos + scale_fill_gradient(low = "blue", high = "red",limits = c(0,max_color))
      }
      p_pos <- p_pos + ggtitle("Positive") +theme(text = element_text(size=9), plot.margin = unit(c(5,1,1,2), "pt"),
                                                  legend.position="bottom",plot.title = element_text(hjust = 0.5)) + theme
      legend <- cowplot::get_legend(p_pos)
      p_pos <- draw_plot(p_pos + theme(legend.position="none"), .5, .15, .5, .85)
    }
    p <-ggdraw() + p_neg + p_pos + draw_grob(legend, 0, 0, 1, .15)
  }
  return(p)
}
