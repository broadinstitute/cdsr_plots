#' Set publication theme on ggplot object
#'
#' @param base_size
#' @param base_family
#'
#' @return
#' @export
#'
#' @examples
#' https://rpubs.com/Koundy/71792
theme_Publication <- function(base_size=12, base_family="Helvetica", base_face = 'bold', legend_bottom = TRUE) {
  library(grid)
  library(ggthemes)
  library(ggplot2)
  if (legend_bottom) {
    loc <- 'bottom'
    orientation <- 'horizontal'
  } else {
    loc <- 'right'
    orientation <- 'vertical'
  }
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = base_face,
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = base_face,size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(),
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = loc,
            legend.text = element_text(size = rel(1.2)),
            legend.direction = orientation,
            legend.key.size= unit(0.3, "cm"),
            legend.spacing = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(5,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face=base_face)
    ))
  
}

#' Set publication fill scheme
#'
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
}

#' Set publication color scheme
#'
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
scale_color_Publication <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
}

