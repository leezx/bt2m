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
FindBifurcationResolution.old <- function(seuratObj, cAlgorithm = 4, resolution.sets = 51, 
                                      graph.name = "RNA_nn", verbose = F) {
  # check
  if (is.null(tmp.seuset.flt@graphs[[graph.name]])) {
    message(sprintf("Input graph: %s is empty, please creat it before...", graph.name))
  }
  # 二分查找
  optimal_resolution <- 0
  # 顺序or倒序
  for (resolution in seq(0,1,length.out = resolution.sets)) {
    # 从小到大搜索到第一个cluster_number>1的resolution，0到多的一个区间
    seuratObj <- FindClusters(seuratObj, algorithm = cAlgorithm, graph.name = graph.name, 
                              resolution = resolution, verbose = F)
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
      if (resolution==0) stop("more than 2 clusters at resolution: 0!!! Please check!!!")
      linearseq <- seq(resolution, resolution-1/(resolution.sets-1), length.out = resolution.sets)
      # 从大到小搜寻第一个等于2的
      for (j in linearseq) {
        seuratObj <- FindClusters(seuratObj, algorithm = cAlgorithm, graph.name = graph.name, 
                                  resolution = j, verbose = F)
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

#' Find the best resolution for bifurcation in graph-based clustering
#'
#' @param seuratObj A Seurat object
#' @param resolution.sets The number of resolution for searching
#'
#' @return The best resolution which can bifurcate all cells
#' @export
#'
FindBifurcationResolution <- function(seuratObj, cAlgorithm = 4, graph.name = "RNA_nn", verbose = F) {
    # check
    if (is.null(tmp.seuset.flt@graphs[[graph.name]])) {
    message(sprintf("Input graph: %s is empty, please creat it before...", graph.name))
    }
    # 二分查找
    max.iter.num <- 50
    tmp.iter.num <- 0
    tmp.min <- 0
    tmp.max <- 1
    tmp.mid <- (tmp.min + tmp.max) / 2
    seuratObj <- FindClusters(seuratObj, algorithm = cAlgorithm, graph.name = graph.name, 
                                  resolution = tmp.mid, verbose = F)
    tmp.cluster.num <- length(unique(seuratObj@active.ident))
    while (tmp.cluster.num != 2) {
        if (tmp.cluster.num < 2) {
            tmp.min <- tmp.mid
            tmp.mid <- (tmp.min + tmp.max) / 2
        } else if (tmp.cluster.num > 2) {
            tmp.max <- tmp.mid
            tmp.mid <- (tmp.min + tmp.max) / 2
        }
        seuratObj <- FindClusters(seuratObj, algorithm = cAlgorithm, graph.name = graph.name, 
                                  resolution = tmp.mid, verbose = F)
        tmp.cluster.num <- length(unique(seuratObj@active.ident))
        tmp.iter.num <- tmp.iter.num + 1
        if (verbose) message(sprintf("mid resolution is: %s, get new cluster number: %s", tmp.mid, tmp.cluster.num))
        if (tmp.iter.num > max.iter.num) {
            stop("fail to find the resolution for 2 clusters")
        }
    }
    return(tmp.mid)
}

#' Bifurcation based on graph-based clustering
#'
#' @param seuratObj A Seurat object
#' @param resolution.sets The number of resolution for searching
#'
#' @return A bifurcated Seurat object (see active.ident).
#' @export
#'
Bt2mBifucation.graph <- function(seuratObj, cAlgorithm = 4, graph.name = "RNA_snn",
                                  slot = "data", assay = "RNA", force = T, verbose = F) {
  # run PCA and SNN
  #message("run Bt2mBifucation.graph...")
  if (verbose) message(paste("Processing ", nrow(seuratObj)," gene and ", ncol(seuratObj), " cells", sep = ""))
  if (force) {
      # redo the PCA and graph, for subgraph
      seuratObj <- RunPCA(seuratObj, slot = slot, assay = assay, verbose = F)
      seuratObj <- FindNeighbors(seuratObj, dims = 1:10, verbose = F)
  }
  # get optimal resolution
  optimal_resolution <- FindBifurcationResolution(seuratObj, cAlgorithm = cAlgorithm, graph.name = graph.name, 
                                                  verbose = verbose)
  # if optimal_resolution==0, indivisiable
  if (optimal_resolution==0) {
    seuratObj@active.ident <- rep(0, ncol(seuratObj))
    names(seuratObj@active.ident) <- colnames(seuratObj)
    return(seuratObj)
  }
  seuratObj <- FindClusters(seuratObj, algorithm = cAlgorithm, graph.name = graph.name, 
                            resolution = optimal_resolution, verbose = F)
  # make sure output clusters are c(0, 1)
  seuratObj@active.ident <- plyr::mapvalues(seuratObj@active.ident, from = unique(seuratObj@active.ident), to = c(0,1))
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
    distMat <- HiClimR::fastCor(t(data.use))
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
