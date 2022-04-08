#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Find the best resolution for bifurcation in graph-based clustering
#' 
#' @param seuratObj A Seurat object
#' @param res_sets The number of resolution for searching
#' 
#' @return The best resolution which can bifurcate all cells
#' @export
#' 
FindBifurcationResolution <- function(seuratObj, res_sets = 50) {
  # 二分查找
  optimal_resolution <- 0
  # 顺序or倒序
  for (resolution in seq(0,1,length.out = res_sets)) {
    # 从小到大搜索到第一个cluster_number>1的resolution，0到多的一个区间
    seuratObj <- FindClusters(seuratObj, resolution = resolution, verbose = F)
    cluster_number <- length(unique(seuratObj@active.ident))
    if (cluster_number>1) {
      # 倒序搜索到第一个
      # must use log order for graph method to do bifurcation, 0.1 may have 3 clusters
      # log divide the interval to make sure we can find the bifurcated resolution parameter
      # for (j in seq(resolution-0.1, resolution, length.out = res_sets)) {
      logseq <- c(1/(res_sets-1) %o% 2^(0:res_sets)) + resolution - 1/(res_sets-1)
      # 从大到小搜寻第一个等于2的
      for (j in logseq) {
        seuratObj <- FindClusters(seuratObj, resolution = j, verbose = F)
        cluster_number <- length(unique(seuratObj@active.ident))
        if (cluster_number==2) {
          optimal_resolution <- j
          break
        }
      }
      break
    }
  }
  return(optimal_resolution)
}

#' Bifurcation based on graph-based clustering
#' 
#' @param seuratObj A Seurat object
#' @param res_sets The number of resolution for searching
#' 
#' @return A bifurcated Seurat object (see active.ident).
#' @export
#' 
IterbiBifucation.graph <- function(seuratObj, res_sets = 50) {
  # run PCA and SNN
  #message("run IterbiBifucation.graph...")
  #message(paste("Processing ", nrow(seuratObj)," gene and ", ncol(seuratObj), " cells", sep = ""))
  seuratObj <- RunPCA(seuratObj, verbose = F)
  seuratObj <- FindNeighbors(seuratObj, dims = 1:10, verbose = F)
  # get optimal resolution
  optimal_resolution <- FindBifurcationResolution(seuratObj, res_sets = res_sets)
  seuratObj <- FindClusters(seuratObj, resolution = optimal_resolution, verbose = F)
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
IterbiBifucation.hclust <- function(seuratObj, method="euclidean") {
  seuratObj <- RunPCA(seuratObj, verbose = F)
  data.use <- Embeddings(object = seuratObj[["pca"]])
  if (method=="euclidean") {
    distMat <- parDist(data.use, threads = 3, method = "euclidean")
  } else if (method=="correlation") {
    distMat <- fastCor(t(data.use))
    distMat <- as.dist(distMat)
  }
  hc <- hclust(distMat)
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
IterbiBifucation.kmeans <- function(seuratObj, method="euclidean") {
  seuratObj <- RunPCA(seuratObj, verbose = F)
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
