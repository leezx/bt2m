#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Prepare expression matrix for heatmap visualization, estimate the overall expression of marker modules
#'
#' @param seuratObj A Seurat object
#' @param bt2m.marker.chain bt2m.marker.chain from bt2m
#' @param assay Assay used for prediction
#' @param slot slot used for prediction
#' @param known_markers the genes to include in the heatmap (like known markers)
#'
#' @return A expression matrix contain the average expression of marker modules and known markers
#' @export
#'
PrepareExpressionMatrix <- function(seuratObj, bt2m.marker.chain, assay="RNA", slot = "scale.data", min.marker.num = 10,
                                              known_markers = c(), verbose = F) {
  # get data
  reference.data <- GetAssayData(object = seuratObj, assay = assay, slot = slot)
  reference.data <- as.matrix(reference.data) # RNA data is not matrix
  if (is.null(reference.data)) stop("Data not exist. Please input right assay and slot!!!")
  all.exprMat <- data.frame()
  tmp.name.list <- c()
  for (i in unique(bt2m.marker.chain$cluster)) {
    tmp.module <- subset(bt2m.marker.chain, cluster==i)$gene
    if (length(tmp.module) < min.marker.num) {
        next
    } else {
        tmp.name.list <- c(tmp.name.list, i)
    }
    if (verbose) message(sprintf("Find %s genes in %s cluster.", length(tmp.module), i))
    tmp.exprMat <- colMeans(reference.data[tmp.module,])
    all.exprMat <- rbind(all.exprMat, tmp.exprMat)
    # break
  }
  rownames(all.exprMat) <- tmp.name.list # unique(bt2m.marker.chain$cluster)
  colnames(all.exprMat) <- colnames(reference.data)
  # gene for validation
  if (length(known_markers)>1) {
    all.exprMat <- rbind(all.exprMat, reference.data[known_markers,])
  }
  return(all.exprMat)
}

#' Colors for bt2m
#'
#' @return A color vector
#' @export
#'
Bt2mColors <- function() {
  bt2m.colors <- unique(c(RColorBrewer::brewer.pal(n = 9, name = "Set1"),
                            RColorBrewer::brewer.pal(n = 8, name = "Dark2"),
                            RColorBrewer::brewer.pal(n = 12, name = "Paired"),
                            RColorBrewer::brewer.pal(n = 8, name = "Set2"),
                            RColorBrewer::brewer.pal(n = 12, name = "Set3"),
                            RColorBrewer::brewer.pal(n = 8, name = "Accent"),
                            RColorBrewer::brewer.pal(n = 9, name = "Pastel1"),
                            RColorBrewer::brewer.pal(n = 8, name = "Pastel2")
  ))
  return(bt2m.colors)
}

#' Draw heatmap by ComplexHeatmap
#'
#' @param seuratObj A Seurat object
#' @param bt2m.cellMeta bt2m.cellMeta from bt2m
#' @param bt2m.marker.chain bt2m.marker.chain from bt2m
#' @param compare_anno select a annotation in seurat object for comparision (e.g. previous annotation)
#' @param known_markers the genes to include in the heatmap (like known markers)
#'
#' @return ComplexHeatmap object
#' @export
#'
DrawMarkerChainHeatmap <- function(seuratObj, bt2m.cellMeta, bt2m.marker.chain, compare_anno="",
                                        known_markers=c(), verbose = F) {
  # get data
  all.exprMat <- PrepareExpressionMatrix(seuratObj, bt2m.marker.chain, slot = "scale.data", 
                                         known_markers=known_markers)
  # heatmap color set
  col_fun = circlize::colorRamp2(c(-1, 0, 1), c("green", "white", "red"))
  #col_fun = circlize::colorRamp2(c(-2, 0, 2), c("#FF00FF","#000000","#FFFF00")) # seurat color
  #col_fun = circlize::colorRamp2(c(-2, 0, 2), c("blue","yellow","red"))
  col_fun(seq(-3, 3))
  # col_fun = colorRampPalette(RColorBrewer::brewer.pal(8, "Blues"))(25) # bad
  # col_fun = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(100)
  # get colors
  bt2m.colors <- Bt2mColors()
  # prepare color
  color.list <- list()
  for (i in paste("L",0:(ncol(bt2m.cellMeta)-1),sep = "")) {
    # print(i)
    tmp.names <- unique(bt2m.cellMeta[[i]])
    tmp.colors <- bt2m.colors[1:length(tmp.names)]
    names(tmp.colors) <- tmp.names
    color.list[[i]] <- tmp.colors
    # break
  }
  if (compare_anno!="") {
    # add our annotation
    tmp.anno <- seuratObj[[compare_anno]]
    tmp.anno.vector <- tmp.anno[,1]
    names(tmp.anno.vector) <- rownames(tmp.anno)
    # add color
    color.list$Annotation <- bt2m.colors[1:length(unique(tmp.anno.vector))]
    names(color.list$Annotation) <- as.character(unique(tmp.anno.vector))
  }
  # sort clusters, to avoid 10 follow 1 cluster
  for (i in colnames(bt2m.cellMeta)) {
    bt2m.cellMeta[,i] <- factor(bt2m.cellMeta[,i], levels = unique(bt2m.cellMeta[,i]))
  }
  # top color bar
  ha <- ComplexHeatmap::HeatmapAnnotation(
    df = bt2m.cellMeta,
    Annotation = tmp.anno.vector[rownames(bt2m.cellMeta)],
    col = color.list
  )
  #
  # options(repr.plot.width=10, repr.plot.height=7)
  if (verbose) message("yes")
  rownames(all.exprMat) <- paste(names(table(bt2m.marker.chain$cluster)[rownames(all.exprMat)]), 
                                 " (",table(bt2m.marker.chain$cluster)[rownames(all.exprMat)]," g)", sep = "")
  ht <- ComplexHeatmap::Heatmap(all.exprMat[,rownames(bt2m.cellMeta)], col = col_fun, name = "Expression",
                #row_split = sc3_marker$cluster, column_split = anno$col$Cluster, layer_fun = layer_fun,
                use_raster = F,
                top_annotation = ha, # left_annotation = la,
                show_row_names = T, show_column_names = F, cluster_rows = F, cluster_columns = F)
  return(ht)
}

#' Draw cluster tree/chain by clustree
#'
#' @param seuratObj A Seurat object
#' @param bt2m.cellMeta bt2m.cellMeta from bt2m
#' @param node_text_size see node_text_size in clustree function
#' @param node_size see node_size in clustree function
#'
#' @return clustree object
#' @export
#'
DrawBt2mClusterTree <- function(seuratObj, bt2m.cellMeta, node_text_size=7, node_size=8) {
  seuratObj <- WriteBt2mIntoSeurat(seuratObj, bt2m.cellMeta)
  #
  # options(repr.plot.width=10, repr.plot.height=10)
  # tree <- clustree(seuratObj, prefix = "L", node_text_size = 7)
  # options(repr.plot.width=4, repr.plot.height=7)
  # options(repr.plot.width=4, repr.plot.height=7)
  tree <- clustree::clustree(seuratObj, prefix = "L", node_text_size = node_text_size, node_size = node_size) +
    theme(legend.position = "none") +
    coord_flip() + scale_y_reverse() +
    scale_color_manual(values=RColorBrewer::brewer.pal(12,"Set3"))
  return(tree)
}

#' Draw cluster tree/chain by clustree
#'
#' @param seuratObj A Seurat object
#' @param bt2m.marker.chain bt2m.marker.chain from bt2m
#' @param rmDup remove duplicated genes or not
#' @param top_n top n genes for dotplot (sort by P-value)
#' @param rel_heights relative heights for plot_grid(), the first and last control the height of top and bottom margin, the middle one controls the height of tree
#' @param rel_widths relative widths for plot_grid(), the first one controls the widths of tree, the second one controls the widths of dotplot
#'
#' @return A plot_grid integrated tree and dotplot
#' @export
#'
DrawMarkerChainDotplot <- function(seuratObj, bt2m.marker.chain, rmDup=T, top_n=3,
                                      rel_heights=c(0.5, 9, 1.5), rel_widths = c(1, 4)) {
  if (rmDup) bt2m.marker.chain.rmDup <- RemoveDuplicatedMarker(bt2m.marker.chain)
  # check level, set to the last
  check_level <- colnames(bt2m.cellMeta)[ncol(bt2m.cellMeta)]
  # set active.ident
  seuratObj@active.ident <- DataframeToVector(seuratObj[[check_level]])
  #
  top <- subset(bt2m.marker.chain.rmDup) %>%
    group_by(cluster) %>%
    dplyr::top_n(top_n, abs(correlation))
  p1 <- DrawBt2mClusterTree(seuratObj, bt2m.cellMeta)
  p2 <- DotPlot(seuratObj, features = rev(top$gene), dot.scale = 7) + RotatedAxis() + ylab(NULL) + xlab(NULL)
  left.anno <- cowplot::plot_grid(NULL,p1,NULL,ncol = 1, rel_heights = rel_heights)
  # options(repr.plot.width=15, repr.plot.height=7)
  p <- cowplot::plot_grid(left.anno, p2, nrow = 1, rel_widths = rel_widths)
  return(p)
}

#' Draw barplot for GO annotation
#'
#' @param barplot_df A GO annotation dataframe from clusterProfiler
#' @param cluster.chain cluster.chain from bt2m
#'
#' @return A ggplot barplot
#' @export
#'
DrawBarplotGO <- function(barplot_df, cluster.chain) {
  # library(Hmisc)
  # library(stringr)
  # library(RColorBrewer)
  colors <- brewer.pal(12, "Set1")
  # colors <- brewer.pal(10,"Paired")
  for (i in 1:dim(barplot_df)[1]) {
    barplot_df[i,]$Description <- capitalize(as.character(barplot_df[i,]$Description))
  }
  # if (length(barplot_df)<2 | is.null(barplot_df)) {next}
  barplot_df$cluster <- factor(barplot_df$cluster, levels=rev(cluster.chain))
  barplot_df <- barplot_df[order(barplot_df$pvalue, decreasing = F),]
  barplot_df <- barplot_df[order(barplot_df$cluster, decreasing = F),]
  maxpvalue <- max(-log10(barplot_df$pvalue))
  # if (length(barplot_df$Description)>20) {barplot_df <- barplot_df[1:20,]}
  g <- ggplot(data=barplot_df, aes(x=Description, y=-log10(pvalue))) +
    # facet_wrap(.~cluster) +
    geom_bar(stat="identity", aes(color=cluster, fill=cluster), alpha=0.8) +
    geom_text(aes(label=Count),color="black",vjust=0.4,hjust=-0.5,size=3,fontface="bold") +
    ylim(0, maxpvalue*1.1) +
    coord_flip() +
    labs(x = "", y = "-Log10(P-value)", title="") +
    theme_bw() +
    # theme(legend.position = "none") +
    theme(axis.text.y = element_text(size = 14, color = "black", face = "plain"),
          axis.text.x = element_text(size = 12, color = "black", face = "plain"),
          axis.title =element_text(size = 15)) +
    # axis.title = element_blank() ,plot.title = element_blank(), axis.ticks.y = element_blank(),
    # axis.ticks.x = element_blank(), axis.text.x = element_blank(),axis.line = element_blank(),
    # panel.border = element_blank(),
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  plot.margin=unit(c(0,0,0,0), "cm"),
          panel.border = element_blank()) +
    theme(axis.line = element_line(color = 'black')) +
    scale_x_discrete(labels=function(x) stringr::str_wrap(x, width=25), limits = rev(barplot_df$Description)) +
    scale_fill_manual(values=colors) +
    scale_color_manual(values=colors)
  g
}

#' Get GO chain corresponding to a cluster chain
#'
#' @param bt2m.GO.anno A GO annotation dataframe from clusterProfiler
#' @param cluster.chain cluster.chain from bt2m
#' @param top_n top n GO terms for barplot (sort by P-value)
#'
#' @return A ggplot barplot showed GO chain
#' @export
#'
DrawGOchain <- function(bt2m.GO.anno, cluster.chain, top_n=3) {
  bt2m.GO.anno <- bt2m.GO.anno[!duplicated(bt2m.GO.anno$ID, fromLast = T),]
  GO.chain <- subset(bt2m.GO.anno, cluster %in% names(cluster.chain))
  top <- subset(GO.chain) %>%
    group_by(cluster) %>%
    top_n(top_n, -qvalue)
  p <- DrawBarplotGO(top, names(cluster.chain))
  return(p)
}

