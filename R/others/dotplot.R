
#' Draw dotplot for human SC3 object
#' @param sce SC3 object
#' @param genes.plot genes
#' @param use group used in SC3 object
#' @param xAngle angle of x title
#' @param cols.use colors used
#' @param col.min col.min
#' @param col.max col.max
#' @param dot.min dot.min
#' @param dot.scale whether to scale or not
#' @param group.order group.order
#' @param scale.by scale.by
#' @param scale.min scale.min
#' @param scale.max scale.max
#' @param group.by group.by
#' @param title title
#' @param plot.legend plot legend or not
#' @param do.return return object or not
#' @param x.lab.rot rotate x title or not
#'
#' @return a ggplot dotplot
#' @export
#' @examples
#' options(repr.plot.width=4, repr.plot.height=6)
#' plot.dotplot.human(sce = sce_HSCR_pure, genes.plot = uniquegenes2, use="cellGroup", group.order = group.order,
#'           scale.min=0, scale.max=100, title="", plot.legend = T, xAngle=90)
#'
plot.dotplot.SC3 <- function (sce, genes.plot, use="cellGroup2",xAngle=60, cols.use = c("lightgrey", "blue"),
                              col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 6, group.order = c(),
                              scale.by = "radius", scale.min = NA, scale.max = NA, group.by, title="",
                              plot.legend = FALSE, do.return = FALSE, x.lab.rot = FALSE)
{
  #cols.use = c("lightgrey", "blue"); genes.plot = tmp$human; sce <- sce_comp; use="cellGroup2";
  #col.min = -2.5; col.max = 2.5; dot.min = 0; dot.scale = 6; group.order = c();xAngle=60;
  #scale.by = "radius"; scale.min = 0; scale.max = 100; title="";
  #plot.legend = T; do.return = FALSE; x.lab.rot = FALSE
  #
  library(dplyr)
  library(tidyr)
  PercentAbove <- function(x, threshold){
    return(length(x = x[x > threshold]) / length(x = x))
  }
  scale.func <- switch(EXPR = scale.by, size = scale_size,
                       radius = scale_radius, stop("'scale.by' must be either 'size' or 'radius'"))
  #if (!missing(x = group.by)) {
  #    object <- SetAllIdent(object = object, id = group.by)
  #}
  #data.to.plot <- data.frame(FetchData(object = object, vars.all = genes.plot))
  data.to.plot <- as.data.frame(t(logcounts(sce)[genes.plot,]))
  colnames(x = data.to.plot) <- genes.plot
  data.to.plot$cell <- rownames(x = data.to.plot)
  #data.to.plot$id <- object@ident
  data.to.plot$id <- colData(sce)[,use]
  data.to.plot <- data.to.plot %>% gather(key = genes.plot,
                                          value = expression, -c(cell, id))
  data.to.plot <- data.to.plot %>% group_by(id, genes.plot) %>%
    summarize(avg.exp = mean(expm1(x = expression)), pct.exp = PercentAbove(x = expression,
                                                                            threshold = 0))
  data.to.plot <- data.to.plot %>% ungroup() %>% group_by(genes.plot) %>%
    dplyr::mutate(avg.exp.scale = scale(x = avg.exp)) %>% dplyr::mutate(avg.exp.scale = MinMax(data = avg.exp.scale,
                                                                                               max = col.max, min = col.min))
  data.to.plot$genes.plot <- factor(x = data.to.plot$genes.plot,
                                    levels = rev(x = genes.plot))
  data.to.plot$pct.exp[data.to.plot$pct.exp < dot.min] <- NA
  data.to.plot$pct.exp <- data.to.plot$pct.exp * 100
  # group.order
  if (length(group.order)>0) {
    data.to.plot$id <- factor(data.to.plot$id, levels = group.order)
  }
  # plot
  p <- ggplot(data = data.to.plot, mapping = aes(y = genes.plot, x = id)) +
    geom_point(mapping = aes(size = pct.exp, color = avg.exp.scale)) +
    scale.func(range = c(0, dot.scale), limits = c(scale.min,scale.max)) +
    theme_bw(base_line_size = 0.1, base_rect_size = 0.1) +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
          panel.border = element_rect(colour = 'black'), panel.grid.minor = element_blank(),
          plot.background = element_blank(), panel.grid.major = element_blank(),
          axis.text.x  = element_text(face="plain", angle=xAngle, size = 10, color = "black", vjust=0.5),
          axis.text.y  = element_text(face="italic", size = 10, color = "black")) +
    labs(title = title) +
    guides(colour = guide_colourbar(title = "Relative expression", title.position = "left",
                                    direction="vertical", title.theme = element_text(angle=90), title.hjust=-4),
           size = guide_legend(title = "Expressed percentage (%)", title.position = "left",
                               direction="vertical", title.theme = element_text(angle=90), title.hjust=-15))
  ####################
  if (length(x = cols.use) == 1) {
    p <- p + scale_color_distiller(palette = cols.use)
  } else {
    p <- p + scale_color_gradient(low = cols.use[1], high = cols.use[2])
  }
  if (!plot.legend) {
    p <- p + theme(legend.position = "none")
  }
  x.lab.rot <- F
  if (x.lab.rot) {
    p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
  }
  suppressWarnings(print(p))
  if (do.return) {
    return(p)
  }
}
