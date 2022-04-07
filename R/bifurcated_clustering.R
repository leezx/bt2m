find_bifurcation_resolution <- function(seuset, res_sets = 50) {
  # 二分查找
  optimal_resolution <- 0
  # 顺序or倒序
  for (resolution in seq(0,1,length.out = res_sets)) {
    # 从小到大搜索到第一个cluster_number>1的resolution，0到多的一个区间
    seuset <- FindClusters(seuset, resolution = resolution, verbose = F)
    cluster_number <- length(unique(seuset@active.ident))
    if (cluster_number>1) {
      # 倒序搜索到第一个
      # must use log order for graph method to do bifurcation, 0.1 may have 3 clusters
      # log divide the interval to make sure we can find the bifurcated resolution parameter
      # for (j in seq(resolution-0.1, resolution, length.out = res_sets)) {
      logseq <- c(1/(res_sets-1) %o% 2^(0:res_sets)) + resolution - 1/(res_sets-1)
      # 从大到小搜寻第一个等于2的
      for (j in logseq) {
        seuset <- FindClusters(seuset, resolution = j, verbose = F)
        cluster_number <- length(unique(seuset@active.ident))
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

iterbi.bifucation.graph <- function(tmp.seuset, res_sets = 50) {
  # run PCA and SNN
  #message("run iterbi.bifucation.graph...")
  #message(paste("Processing ", nrow(tmp.seuset)," gene and ", ncol(tmp.seuset), " cells", sep = ""))
  tmp.seuset <- RunPCA(tmp.seuset, verbose = F)
  tmp.seuset <- FindNeighbors(tmp.seuset, dims = 1:10, verbose = F)
  # get optimal resolution
  optimal_resolution <- find_bifurcation_resolution(tmp.seuset, res_sets = res_sets)
  tmp.seuset <- FindClusters(tmp.seuset, resolution = optimal_resolution, verbose = F)
  return(tmp.seuset)
}

iterbi.bifucation.hclust <- function(tmp.seuset, method="euclidean") {
  tmp.seuset <- RunPCA(tmp.seuset, verbose = F)
  data.use <- Embeddings(object = tmp.seuset[["pca"]])
  if (method=="euclidean") {
    distMat <- parDist(data.use, threads = 3, method = "euclidean")
  } else if (method=="correlation") {
    distMat <- fastCor(t(data.use))
    distMat <- as.dist(distMat)
  }
  hc <- hclust(distMat)
  binaryC <- cutree(hc, k = 2)
  tmp.seuset@active.ident <- factor(binaryC[names(tmp.seuset@active.ident)]-1)
  return(tmp.seuset)
}

iterbi.bifucation.kmeans <- function(tmp.seuset, method="euclidean") {
  tmp.seuset <- RunPCA(tmp.seuset, verbose = F)
  data.use <- Embeddings(object = tmp.seuset[["pca"]])
  if (method=="euclidean") {
    distMat <- parDist(data.use, threads = 3, method = "euclidean")
  } else if (method=="correlation") {
    distMat <- fastCor(t(data.use))
    distMat <- as.dist(distMat)
  }
  #hc <- hclust(distMat)
  #binaryC <- cutree(hc, k = 2)
  binaryC <- kmeans(distMat, 2)$cluster
  tmp.seuset@active.ident <- factor(binaryC[names(tmp.seuset@active.ident)]-1)
  return(tmp.seuset)
}
