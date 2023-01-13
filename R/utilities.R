
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Add mising genes to the gene expression matrix (for multiple datasets integraiton)
#'
#' @param epxrM expression matrix in data.frame, matrix or dgCMatrix format
#' @param all.genes full list of genes
#'
#' @return A gene expression matrix in dgCMatrix format
#' @export
#'
add.missing.genes <- function(epxrM, all.genes) {
  # to matrix
  epxrM <- as.matrix(epxrM)
  epxrM <- epxrM[rownames(epxrM) %in% all.genes,]
  # # check how many genes are not detected
  missing.genes <- all.genes[!all.genes %in% rownames(epxrM)]
  # add missing expr to raw count data is good for futher integration
  add.df <- matrix(rep(0, times = length(missing.genes), each = ncol(epxrM)),
                   nrow = length(missing.genes),
                   ncol = ncol(epxrM))
  rownames(add.df) <- missing.genes
  colnames(add.df) <- colnames(epxrM)
  epxrM <- rbind(epxrM, add.df)
  epxrM <- epxrM[all.genes,]
  epxrM <- as(epxrM, "dgCMatrix")
  return(epxrM)
}

# facet with background color
# path: project/Data_center/analysis/ApcKO_multiomics/ApcKO-seurat.ipynb
# add full cells as background point to show skeleton
add.background.point.facet.wrap <- function(all_tsne, cluster.order, sample.order,
                                       cluster="cluster", sample="sample", dim1="UMAP_1", dim2="UMAP_2"
) {
  # project/Data_center/analysis/ApcKO_multiomics/ApcKO-seurat.ipynb
  # add full cells as background point to show skeleton
  # repeat the col or row for n times
  rep.row<-function(x,n){
    matrix(rep(x,each=n),nrow=n)
  }
  rep.col<-function(x,n){
    matrix(rep(x,each=n), ncol=n, byrow=TRUE)
  }
  # bind_rows(replicate(3, all_tsne, simplify = F))
  samples <- unique(as.character(all_tsne[,sample]))
  # replicate each samples in coloum
  rep.df <- rep.col(all_tsne[,cluster], length(samples))
  colnames(rep.df) <- samples
  # add to origin
  all_tsne_rep <- cbind(all_tsne, rep.df)
  # each row is a cell,
  # if the cell is not used in facet, then assign it to others
  for (i in samples) {
    all_tsne_rep[all_tsne[,sample]!=i, i] <- "others"
    next
  }
  all_tsne_rep <- all_tsne_rep[,c(dim1, dim2, samples)]
  # melt to add full list of samples, cell number duplicated for length(samples)
  merged_df <- reshape2::melt(all_tsne_rep, id.vars = c(dim1, dim2))
  # set levels for samples
  merged_df$variable <- factor(merged_df$variable, levels = sample.order)
  # set levels for clusters
  merged_df$value <- factor(merged_df$value, levels = c(cluster.order, "others"))
  merged_df <- merged_df[order(merged_df$value, decreasing = T),]
  colnames(merged_df) <- c(dim1, dim2, sample, cluster)
  return(merged_df)
}

# three variable, use facet_grid
# path: project/Data_center/analysis/ApcKO_multiomics/ApcKO-seurat.ipynb
add.background.point.facet.grid <- function(all_tsne, cluster.order, sample.order, assay.order,
                                            cluster="cluster", sample="sample", assay="assay",
                                            dim1="UMAP_1", dim2="UMAP_2"
) {
  # project/Data_center/analysis/ApcKO_multiomics/ApcKO-seurat.ipynb
  # add full cells as background point to show skeleton
  # repeat the col or row for n times
  rep.row<-function(x,n){
    matrix(rep(x,each=n),nrow=n)
  }
  rep.col<-function(x,n){
    matrix(rep(x,each=n), ncol=n, byrow=TRUE)
  }
  # bind_rows(replicate(3, all_tsne, simplify = F))
  samples <- unique(as.character(all_tsne[,sample]))
  # replicate each samples in coloum
  rep.df <- rep.col(all_tsne[,cluster], length(samples))
  colnames(rep.df) <- samples
  # add to origin
  all_tsne_rep <- cbind(all_tsne, rep.df)
  # each row is a cell,
  # if the cell is not used in facet, then assign it to others
  for (i in samples) {
    all_tsne_rep[all_tsne[,sample]!=i, i] <- "others"
    next
  }

  all_tsne_rep <- all_tsne_rep[,c(dim1, dim2, samples, assay)]
  # melt to add full list of samples, cell number duplicated for length(samples)
  merged_df <- reshape2::melt(all_tsne_rep, id.vars = c(dim1, dim2, assay))
  # set levels for samples
  merged_df$variable <- factor(merged_df$variable, levels = sample.order)
  # set levels for clusters
  merged_df$value <- factor(merged_df$value, levels = c(cluster.order, "others"))
  merged_df <- merged_df[order(merged_df$value, decreasing = T),]
  colnames(merged_df) <- c(dim1, dim2, assay, sample, cluster)
  merged_df$assay <- factor(merged_df$assay, levels = assay.order)
  return(merged_df)
}


#' A general function to identify DEGs between case and control
DEG.cluster.list <- function(seuratObj, cluster.list, ident.1 = "Vcl cKO", ident.2 = "Control", assay = "RNA") {
    DEGs <- list()
    for (i in names(cluster.list)) {
        condition <- seuratObj$cluster %in% cluster.list[[i]]
        tmp.seuratObj <- subset_cells(seuratObj, condition)
        tmp.seuratObj@active.ident <- tmp.seuratObj$group
        tmp.DEGs <- FindMarkers(tmp.seuratObj, ident.1 = ident.1, ident.2 = ident.2, only.pos = F,
                                min.pct = 0, min.diff.pct = "-Inf", logfc.threshold = 0, assay = assay)
        DEGs[[i]] <- add.missing.DEGs(tmp.DEGs, rownames(tmp.seuratObj@assays$RNA@counts))
        # fill emplty genes
        tmp.empty <- rowSums(DEGs[[i]])==0
        DEGs[[i]][tmp.empty,]$p_val <- 1
        DEGs[[i]][tmp.empty,]$p_val_adj <- 1
        # add gene column
        DEGs[[i]]$gene <- rownames(DEGs[[i]])
        # add log2FC and correlation
        cells_1 <- rownames(subset(seuratObj@meta.data, cluster %in% cluster.list[[i]] & group==ident.2))
        cells_2 <- rownames(subset(seuratObj@meta.data, cluster %in% cluster.list[[i]] & group==ident.1))
        log2fc <- log2(rowMeans(seuratObj@assays$RNA@data[,cells_2])+1) - log2(rowMeans(seuratObj@assays$RNA@data[,cells_1])+1)
        corM <- cor(t(as.matrix(seuratObj@assays$RNA@data[,c(cells_1,cells_2)])), c(rep(0, length.out = length(cells_1)),
                                                                       rep(1, length.out = length(cells_2))
                                                                       ))
        corM[is.na(corM)] <- 0
        DEGs[[i]]$log2FC <- log2fc[DEGs[[i]]$gene]
        DEGs[[i]]$correlation <- as.data.frame(corM)[DEGs[[i]]$gene,]
        # break
    }
    return(DEGs)
}



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions for plotting
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





