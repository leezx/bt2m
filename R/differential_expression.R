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
IdentifyBinaryMarkers <- function(seuratObj, p_value=0.01, min_avg_logFC=0.1, min_correlation=0.3){
  # need high quality markers, correlation is the best indicator
  binary_markers <- list()
  # Positive values indicate that the gene is more highly expressed in the first group
  tmp.markers <- FindMarkers(seuratObj, ident.1 = 0, ident.2 = 1, min.pct = 0.25, only.pos = F, verbose = F)
  # reversed the avg_logFC, set group 1 as control (base-line)
  tmp.markers$avg_logFC <- -tmp.markers$avg_logFC
  tmp.markers$gene <- rownames(tmp.markers)
  # add correlation
  pseudo.phenotype <- plyr::mapvalues(as.integer(seuratObj@active.ident), from = 1:2, to = c(-1,1))
  corMat <- cor(t(as.matrix(seuratObj@assays$RNA@data)), pseudo.phenotype)
  tmp.markers$correlation <- corMat[tmp.markers$gene,]
  #
  tmp.markers <- subset(tmp.markers, p_val<p_value & p_val_adj<p_value & (avg_logFC*correlation > 0))
  # strong filtering
  tmp.markers <- subset(tmp.markers, abs(avg_logFC)>min_avg_logFC, abs(correlation)>min_correlation)
  # b1 is the marker for group 1
  tmp.b1.markers <- subset(tmp.markers, avg_logFC<0)
  # b2 is the marker for group 2
  tmp.b2.markers <- subset(tmp.markers, avg_logFC>0)
  binary_markers[["b1"]] <- tmp.b1.markers
  binary_markers[["b2"]] <- tmp.b2.markers
  return(binary_markers)
}

