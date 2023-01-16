#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Identify binary markers
#' 
#' @param seuratObj A Seurat object
#' @param p_value The number of resolution for searching
#' @param min_avg_logFC Minimal threshold for logFC to pick markers
#' @param min_correlation Minimal threshold for correlation to pick markers
#' 
#' @return A list. b1 is the positive marker of cluster 0, b2 is the positive marker of cluster 1.
#' @export
#' 
FindBinaryMarkers <- function(seuratObj, method = "presto", inputFactor = F, slot = "data", assay = "RNA", 
                              p_value = 0.05, min_avg_log2FC = 0.1, min_correlation = 0.2, verbose = F){
  # first check
  if (length(unique(seuratObj@active.ident)) != 2) stop("The cluster number is not 2 for active.ident!!!")
  # check if the input clusters are factor or not
  if (inputFactor) {
    tmp.levels <- levels(seuratObj@active.ident)
    if (is.null(tmp.levels)) stop("Please transform the cluster name to factors!!!")
    if (length(tmp.levels) != 2) stop("Input factors, but the level number is not 2 for active.ident!!!")
    seuratObj@active.ident <- plyr::mapvalues(seuratObj@active.ident, from = tmp.levels, to = c(0,1))
  } else {
    if (sum(sort(unique(seuratObj@active.ident)) == c(0,1)) != 2) stop("The cluster names are not 0 and 1.")
  }
  # need high quality markers, correlation is the best indicator
  binary_markers <- list()
  # Positive values indicate that the gene is more highly expressed in the first group
  if (method == "presto") {
    # very fast
    seuratObj$seurat_clusters <- seuratObj@active.ident
    if (verbose) {message(sprintf("FindBinaryMarkers for cluster: %s.", paste(unique(seuratObj$seurat_clusters, collapse=", "))))}
    tmp.markers <- presto:::wilcoxauc.Seurat(seuratObj, group_by = 'seurat_clusters', 
                                             assay = slot, seurat_assay = assay)
    tmp.markers <- subset(tmp.markers, group==0)
    rownames(tmp.markers) <- tmp.markers$feature
    tmp.markers <- tmp.markers[,c("pval", "logFC", "pct_in", "pct_out", "padj")]
    colnames(tmp.markers) <- c('p_val','avg_log2FC','pct_in','pct_out','p_val_adj')
    tmp.markers$pct_in <- tmp.markers$pct_in/100
    tmp.markers$pct_out <- tmp.markers$pct_out/100
    tmp.markers <- tmp.markers[order(tmp.markers$p_val_adj, decreasing = F),]
  } else if (method == "seurat") {
      # too slow
      tmp.markers <- FindMarkers(seuratObj, ident.1 = 0, ident.2 = 1, min.pct = 0.25, only.pos = F, verbose = F)
      colnames(tmp.markers) <- c('p_val','avg_log2FC','pct_in','pct_out','p_val_adj')
  }
  if (verbose) message(sprintf("The number of raw marker is %s.", nrow(tmp.markers)))
  # reversed the avg_log2FC, set group 1 as control (base-line)
  tmp.markers$avg_log2FC <- -tmp.markers$avg_log2FC
  tmp.markers$gene <- rownames(tmp.markers)
  # tmp.markers$`pct.1` <- NULL
  # tmp.markers$`pct.2` <- NULL
  # add correlation
  pseudo.phenotype <- plyr::mapvalues(as.integer(seuratObj@active.ident), from = 1:2, to = c(-1,1))
  corMat <- cor(t(as.matrix(seuratObj@assays$RNA@data)), pseudo.phenotype)
  tmp.markers$correlation <- corMat[tmp.markers$gene,]
  #
  tmp.markers <- subset(tmp.markers, p_val<p_value & p_val_adj<p_value & (avg_log2FC*correlation > 0))
  # strong filtering
  tmp.markers <- subset(tmp.markers, abs(avg_log2FC)>min_avg_log2FC & abs(correlation)>min_correlation)
  # b1 is the marker for group 1
  tmp.b1.markers <- subset(tmp.markers, avg_log2FC<0)
  # reverse negative to positive value
  tmp.b1.markers$avg_log2FC <- -tmp.b1.markers$avg_log2FC
  tmp.b1.markers$correlation <- -tmp.b1.markers$correlation
  # b2 is the marker for group 2
  tmp.b2.markers <- subset(tmp.markers, avg_log2FC>0)
  #
  binary_markers[["b1"]] <- tmp.b1.markers
  binary_markers[["b2"]] <- tmp.b2.markers
  return(binary_markers)
}

