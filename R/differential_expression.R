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
FindBinaryMarkers <- function(seuratObj, method = "presto", slot = "data", assay = "RNA", 
                              p_value = 0.05, min_avg_log2FC = 0.1, min_correlation = 0.2){
  # need high quality markers, correlation is the best indicator
  binary_markers <- list()
  # Positive values indicate that the gene is more highly expressed in the first group
  if (method == "presto") {
    # very fast
    tmp.seuset.flt$seurat_clusters <- tmp.seuset.flt@active.ident
    tmp.markers <- presto:::wilcoxauc.Seurat(tmp.seuset.flt, group_by = 'seurat_clusters', 
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

