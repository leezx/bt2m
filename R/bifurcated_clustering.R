#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Find the best resolution for bifurcation in graph-based clustering
#'
#' @param seuratObj A Seurat object
#' @param resolution.sets The number of resolution for searching
#'
#' @return The best resolution which can bifurcate all cells
#' @export
#'
FindBifurcationResolution <- function(seuratObj, resolution.sets = 51, graph.name = "RNA_snn", verbose = F) {
  # check
  if (is.null(tmp.seuset.flt@graphs[[graph.name]])) {
    message(sprintf("Input graph: %s is empty, please creat it before...", graph.name))
  }
  # 二分查找
  optimal_resolution <- 0
  # 顺序or倒序
  for (resolution in seq(0,1,length.out = resolution.sets)) {
    # 从小到大搜索到第一个cluster_number>1的resolution，0到多的一个区间
    seuratObj <- FindClusters(seuratObj, graph.name = graph.name, resolution = resolution, verbose = F)
    cluster_number <- length(unique(seuratObj@active.ident))
    if (verbose) message(sprintf("resolution: %s -> cluster number: %s", resolution, cluster_number))
    if (cluster_number>1) {
      if (verbose) message(sprintf("cluster_number>1 at resolution: %s", resolution))
      # 倒序搜索到第一个
      # must use log order for graph method to do bifurcation, 0.1 may have 3 clusters
      # log divide the interval to make sure we can find the bifurcated resolution parameter
      # fail for discrete data
      # logseq <- c(1/(resolution.sets-1) %o% 2^(0:resolution.sets)) + resolution - 1/(resolution.sets-1)
      # would be fail for continous data, sometimes close to resolution, sometimes close to 0
      # linearseq <- c(seq(resolution,resolution*4/5,length.out = resolution.sets), seq(resolution*4/5,0,length.out = resolution.sets))
      linearseq <- seq(resolution, resolution-1/(resolution.sets-1), length.out = resolution.sets)
      # 从大到小搜寻第一个等于2的
      for (j in linearseq) {
        seuratObj <- FindClusters(seuratObj, graph.name = graph.name, resolution = j, verbose = F)
        cluster_number <- length(unique(seuratObj@active.ident))
        if (verbose) message(sprintf("resolution: %s -> cluster number: %s", j, cluster_number))
        if (cluster_number==2) {
          if (verbose) message(sprintf("cluster_number==2 at resolution: %s", j))
          optimal_resolution <- j
          return(optimal_resolution)
        }
      }
      break
    }
  }
}

#' Bifurcation based on graph-based clustering
#'
#' @param seuratObj A Seurat object
#' @param resolution.sets The number of resolution for searching
#'
#' @return A bifurcated Seurat object (see active.ident).
#' @export
#'
Bt2mBifucation.graph <- function(seuratObj, graph.name = "RNA_snn", resolution.sets = 50, slot = "data", assay = "RNA", force = T, verbose = F) {
  # run PCA and SNN
  #message("run Bt2mBifucation.graph...")
  if (verbose) message(paste("Processing ", nrow(seuratObj)," gene and ", ncol(seuratObj), " cells", sep = ""))
  if (force) {
      # redo the PCA and graph, for subgraph
      seuratObj <- RunPCA(seuratObj, slot = slot, assay = assay, verbose = F)
      seuratObj <- FindNeighbors(seuratObj, dims = 1:10, verbose = F)
  }
  # get optimal resolution
  optimal_resolution <- FindBifurcationResolution(seuratObj, graph.name = graph.name, resolution.sets = resolution.sets, verbose = verbose)
  # if optimal_resolution==0, indivisiable
  if (optimal_resolution==0) {
    seuratObj@active.ident <- rep(0, ncol(seuratObj))
    names(seuratObj@active.ident) <- colnames(seuratObj)
    return(seuratObj)
  }
  seuratObj <- FindClusters(seuratObj, graph.name = graph.name, resolution = optimal_resolution, verbose = F)
  return(seuratObj)
}

#' Bifurcation based on hierarchical clustering
#'
#' @param seuratObj A Seurat object
#' @param method The distance measurement method "euclidean or correlation"
#'
#' @return A bifurcated Seurat object (see active.ident).
#' @export
#'
Bt2mBifucation.hclust <- function(seuratObj, method="euclidean", slot = "data", assay = "RNA") {
  seuratObj <- RunPCA(seuratObj, slot = slot, assay = assay, verbose = F)
  data.use <- Embeddings(object = seuratObj[["pca"]])
  if (method=="euclidean") {
    distMat <- parDist(data.use, threads = 3, method = "euclidean")
  } else if (method=="correlation") {
    distMat <- fastCor(t(data.use))
    distMat <- as.dist(distMat)
  }
  hc <- fastcluster::hclust(distMat)
  binaryC <- cutree(hc, k = 2)
  seuratObj@active.ident <- factor(binaryC[names(seuratObj@active.ident)]-1)
  return(seuratObj)
}

#' Bifurcation based on K-means clustering
#'
#' @param seuratObj A Seurat object
#' @param method The distance measurement method "euclidean or correlation"
#'
#' @return A bifurcated Seurat object (see active.ident).
#' @export
#'
Bt2mBifucation.kmeans <- function(seuratObj, method="euclidean", slot = "data", assay = "RNA") {
  seuratObj <- RunPCA(seuratObj, slot = slot, assay = assay, verbose = F)
  data.use <- Embeddings(object = seuratObj[["pca"]])
  if (method=="euclidean") {
    distMat <- parDist(data.use, threads = 3, method = "euclidean")
  } else if (method=="correlation") {
    distMat <- fastCor(t(data.use))
    distMat <- as.dist(distMat)
  }
  #hc <- hclust(distMat)
  #binaryC <- cutree(hc, k = 2)
  binaryC <- kmeans(distMat, 2)$cluster
  seuratObj@active.ident <- factor(binaryC[names(seuratObj@active.ident)]-1)
  return(seuratObj)
}
