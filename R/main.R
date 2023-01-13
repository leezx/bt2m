#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
globalVariables(
  names=c("Count","Description","avg_logFC","bcg_pct","cluster","cluster_pct","correlation","bt2m.cellMeta","bt2m.marker.chain","markerList","new_cluster","old_cluster","p_val","p_val_adj","pvalue","qvalue","select_cells","seuratObj"),
  package = "bt2m",
  add = TRUE
)

# new functions
# judge the discrete and continuous cell
discreteOrContinuous <- function(seuratObj, group) {
  #
  return("discrete or continuous")
}

# find best marker along 1 chain
bestMarkerAlongChain <- function(seuratObj, moduleChain) {
  #
  return(moduleChain.rmDup)
}

# judge if a cell group is divisible or not
cellDivisibilityCheck <- function(seuratObj) {
  # no high-quality unique marker
  # too many depth for contious cell group
  return(seuratObj)
}

#' Nitialize feature annotation data for each assay
#'
#' @param seuratObj A Seurat object
#' @param slot The method to perform bifurcation clustering "graph (default), hclust or kmeans"
#' @param assay Minimal number of markers to confirm a bifurcation
#'
#' @return A list. cellMeta contains the preliminary bifurcation for each level
#' marker_chain contains all the significant markers for each cluster
#' bifucation contains the bifurcation details (parent, child1, child2)
#' @export
#'
InitializeFeatureAnno <- function(seuratObj, slot = "count", assay = "RNA") {
    tmp.exprM <- GetAssayData(seuratObj, slot = slot, assay = assay)
    seuratObj@misc$feature.data <- list()
    seuratObj@misc$feature.data[[assay]] <- data.frame(row.names = rownames(seuratObj@assays$RNA@counts),
                                                   gene_name = rownames(seuratObj@assays$RNA@counts))
    seuratObj
}

#' The main function to perform iteratively bifurcation clustering
#'
#' @param seuratObj A Seurat object
#' @param method The method to perform bifurcation clustering "graph (default), hclust or kmeans"
#' @param min.marker.num Minimal number of markers to confirm a bifurcation
#' @param max.level.num Maximum number of level for bifurcation
#' @param min.cell.count Minimal number of cells to perform bifurcation (must bigger than PC number: 50)
#' @param resolution.sets The number of resolution for searching
#' @param verbose Print detail proccessing messages
#'
#' @return A list. cellMeta contains the preliminary bifurcation for each level
#' marker_chain contains all the significant markers for each cluster
#' bifucation contains the bifurcation details (parent, child1, child2)
#' @export
#'
RunBt2m.single <- function(seuratObj, slot = "data", assay = "RNA", method = "graph", min.marker.num = 100, 
                    max.level.num = 20, min.cell.count = 50, resolution.sets = 51, verbose = T) {
  # basic info
  message(paste("Current Default Assay is", DefaultAssay(seuratObj), sep = " "))
  seuratObj$cellName <- colnames(seuratObj)
  # check info
  tmp.exprM <- GetAssayData(seuratObj, slot = slot, assay = assay)
  if (length(tmp.exprM) == 0) {
    message(sprintf("Input data under slot: %s, and assay: %s, is empty, please creat/scale it...", slot, assay))
  } else {
    message(sprintf("Input data under slot: %s, and assay: %s, is %s cell x %s features", 
                    slot, assay, ncol(tmp.exprM), nrow(tmp.exprM)))
  }
  # key index
  # level index: 1-20
  # cluster index: L1_(1..n)
  # method = "graph"; min.marker.num = 100; max.level.num = 20; min.cell.count = 50; verbose = T # for quite test
  # suppress warnings
  options(warn = -1, stringsAsFactors = F)
  # min.marker.num <- 100
  # max.level.num <- 20
  # min.cell.count <- 50 # must be more than number PCs
  # verbose <- T
  # key file 1: the cell annotation dataframe
  seuratObj$cellName <- as.character(colnames(seuratObj))
  # create meta data
  bt2m.cellMeta <- data.frame(row.names = seuratObj$cellName, cellName=seuratObj$cellName, stringsAsFactors = F)
  # add initial one
  bt2m.cellMeta[,"L0"] <- 1
  # add default L1-max.level.num
  for (i in paste("L", 1:max.level.num, sep = "")) {
    # print(i)
    bt2m.cellMeta[,i] <- 0
  }
  #
  # cannot start from specifc level, the marker list will be chaos
  bt2m.marker.chain <- data.frame()
  bt2m.bifucation <- data.frame()
  # iteratively bifurcation until no enough markers can be found
  for (i in 0:max.level.num) {
    # prevoius.level.index <- paste("L",i-1,sep = "")
    tmp.level.index <- paste("L",i,sep = "") # for getting data
    next.level.index <- paste("L",i+1,sep = "") # for writing data
    message(paste("We are now at", tmp.level.index, sep = " "))
    # small testing
    # if (tmp.level.index=="L3") break
    #
    for (j in unique(bt2m.cellMeta[,tmp.level.index])) {
      # print(j)
      tmp.cluster.index <- paste(tmp.level.index, j, sep = "_")
      # get cells
      # message("get tmp.cells...")
      tmp.cells <- bt2m.cellMeta[bt2m.cellMeta[,tmp.level.index]==j,]$cellName
      # message(length(tmp.cells))
      # if there are no marker for the previous level, endpoint
      if (grepl("end", j)) {
        bt2m.cellMeta[tmp.cells, next.level.index] <- j
        # message(paste(tmp.cluster.index, "was an end point, skipping...", sep = " "))
        next
      } else {
        if (verbose) message(paste("Bifurcating", tmp.cluster.index, "...", sep = " "))
      }
      # if one column has "end point", it will be characters, need transform it
      j <- as.numeric(j) # proccess the end point
      if (length(tmp.cells) <= min.cell.count) {
        # exit when cells is less than 50
        message(paste("only have ", length(tmp.cells), " cells, set it as an end node", sep=""))
        bt2m.cellMeta[tmp.cells, next.level.index] <- paste("end", next.level.index, 2*j-1, sep = "-")
        # message(paste(tmp.cluster.index, "has no enough cells, this is a end point...", sep = " "))
        next
      }
      ###################
      # binary clustering
      # message(length(tmp.cells))
      # subset funciton have problem
      tmp.seuratObj <- subsetSeuratObjByCells(seuratObj, tmp.cells)
      # tmp.seuratObj <- subset(seuratObj, subset = cellName %in% tmp.cells) 
      # can't find tmp.cells, don't know why?
      #
      if (method == "graph") {
          tmp.seuratObj <- Bt2mBifucation.graph(tmp.seuratObj, resolution.sets = resolution.sets, 
                                                slot = slot, assay = assay)
      } else if (method == "hclust") {
          tmp.seuratObj <- Bt2mBifucation.hclust(tmp.seuratObj, slot = slot, assay = assay)
      } else if (method == "kmeans") {
          tmp.seuratObj <- Bt2mBifucation.kmeans(tmp.seuratObj, slot = slot, assay = assay)
      } else {message("Please select one method in: graph, hclust, kmeans!")}
      # check clustering result
      if (length(unique(tmp.seuratObj@active.ident)) != 2) {
        message("Cannot do bifurcation. Set bigger resolution.sets!!!")
        stop()
      }
      # assign name to the marker dfs
      next.cluster.index.1 <- paste(next.level.index, 2*j-1, sep = "_")
      next.cluster.index.2 <- paste(next.level.index, 2*j, sep = "_")
      # write annotation to the next level
      bt2m.cellMeta[tmp.seuratObj$cellName, next.level.index] <- as.integer(tmp.seuratObj@active.ident) + 2*(j-1)
      # get markers
      binary_markers <- FindBinaryMarkers(tmp.seuratObj)
      b1.markers <- binary_markers[["b1"]]
      b2.markers <- binary_markers[["b2"]]
      # set end point
      if (nrow(b1.markers) >= min.marker.num) {
        b1.markers$parent <- tmp.cluster.index
        b1.markers$cluster <- next.cluster.index.1
        bt2m.marker.chain <- rbind(bt2m.marker.chain, b1.markers)
      }
      if (nrow(b2.markers) >= min.marker.num) {
        b2.markers$parent <- tmp.cluster.index
        b2.markers$cluster <- next.cluster.index.2
        bt2m.marker.chain <- rbind(bt2m.marker.chain, b2.markers)
      }
      # final judgement
      if ((nrow(b1.markers) < min.marker.num) && (nrow(b2.markers) < min.marker.num)) {
        # message(paste(tmp.cluster.index, "has no marker, this is a end point...", sep = " "))
        # set end point (take care, what the name is, when no marker)
        # must add "end" and "level" to avoid the mixture between levels
        bt2m.cellMeta[tmp.cells, next.level.index] <- paste("end", next.level.index, 2*j-1, sep = "-")
      } else {
        tmp.bifucation <- data.frame(parent=tmp.cluster.index, child1=next.cluster.index.1, child2=next.cluster.index.2)
        bt2m.bifucation <- rbind(bt2m.bifucation, tmp.bifucation)
        if (verbose) message(sprintf("Successfully split %s to %s and %s", tmp.cluster.index,
                                     next.cluster.index.1, next.cluster.index.2))
      }
    }
    # exist running when all cluster are end points
    if (sum(!grepl("end", unique(bt2m.cellMeta[,tmp.level.index]))) == 0) {
      message("Bifurcating stopped! No more clusters can be split")
      break
    }
  }
  # remove emplty levels
  bt2m.cellMeta <- bt2m.cellMeta[,colSums(bt2m.cellMeta!=0)!=0]
  bt2m.cellMeta$cellName <- NULL
  # remove same columns
  bt2m.cellMeta <- bt2m.cellMeta[,!duplicated(t(bt2m.cellMeta))]
  # sort the df last to first, not the final one
  for (i in ncol(bt2m.cellMeta):1) {
    bt2m.cellMeta <- bt2m.cellMeta[order( bt2m.cellMeta[,i] ),]
  }
  return(list(cellMeta=bt2m.cellMeta, marker_chain=bt2m.marker.chain, bifucation=bt2m.bifucation))
}

#' order the clusters according to similarity inside the bt2m result
#'
#' @param seuratObj A Seurat object
#' @param bt2m.result A result file from RunBt2m() function

#' @return A re-ordered list. cellMeta contains the final bifurcation for each level
#' marker_chain contains all the significant markers for each cluster
#' bifucation contains the bifurcation details (parent, child1, child2)
#' @export
#'
OrderCluster <- function(seuratObj, bt2m.result) {
    bt2m.cellMeta <- bt2m.result[["cellMeta"]]
    bt2m.marker.chain <- bt2m.result[["marker_chain"]]
    bt2m.bifucation <- bt2m.result[["bifucation"]]
    #
    all.exprMat <- PrepareExpressionMatrix(seuratObj, bt2m.marker.chain, assay = "RNA", slot = "data")
    all.exprMat <- as.data.frame(t(all.exprMat))
    all.exprMat <- all.exprMat[rownames(bt2m.cellMeta),]
    #
    # too hard, need to check the switch logic, compare the former and latter, it's very complicated
    for (i in 1:nrow(bt2m.bifucation)) {
        # print(i)
        child1 <- bt2m.bifucation[i,2]
        child2 <- bt2m.bifucation[i,3]
        tmp.level <- strsplit(child1, split = "_")[[1]][1]
        tmp.level.index <- as.numeric(substr(tmp.level, 2, 10))
        if (tmp.level.index==1) next
        child1.index <- as.numeric(strsplit(child1, split = "_")[[1]][2])
        child2.index <- as.numeric(strsplit(child2, split = "_")[[1]][2])
        #
        tmp.position <- unique(bt2m.cellMeta[,tmp.level])
        # message(paste(tmp.position, collapse = ", "))
        former.child <- tmp.position[tmp.position %in% c(child1.index, child2.index)][1]
        latter.child <- tmp.position[tmp.position %in% c(child1.index, child2.index)][2]
        former.child.index <- which(tmp.position==former.child)
        latter.child.index <- which(tmp.position==latter.child)
        # distance
        all.exprMat$level <- bt2m.cellMeta[,tmp.level]
        module.exprMat <- aggregate(all.exprMat[, 1:(ncol(all.exprMat)-1)], list(all.exprMat$level), mean)
        rownames(module.exprMat) <- module.exprMat$`Group.1`
        module.exprMat$`Group.1` <- NULL
        disMat <- cor(t(module.exprMat), method = "spearman")
        disMat <- 1 - disMat
        #
        if (former.child.index==1) {
            # no formal clusters
            latter.clusters <- tmp.position[(latter.child.index+1):length(tmp.position)]
            former.min.dist <- min(disMat[as.character(former.child), as.character(latter.clusters)])
            latter.min.dist <- min(disMat[as.character(latter.child), as.character(latter.clusters)])
            if (former.min.dist < latter.min.dist) {
                # switch
                tmp.position[former.child.index] <- latter.child
                tmp.position[latter.child.index] <- former.child
                message(sprintf("switch %s %s at level %s", former.child, latter.child, tmp.level))
            }
        } else if (latter.child.index==length(tmp.position)) {
            # no latter clusters
            former.clusters <- tmp.position[1:(former.child.index-1)]
            former.min.dist <- min(disMat[as.character(former.child), as.character(former.clusters)])
            latter.min.dist <- min(disMat[as.character(latter.child), as.character(former.clusters)])
            if (former.min.dist > latter.min.dist) {
                # switch
                tmp.position[former.child.index] <- latter.child
                tmp.position[latter.child.index] <- former.child
                message(sprintf("switch %s %s at level %s", former.child, latter.child, tmp.level))
            }
        } else {
            former.clusters <- tmp.position[1:(former.child.index-1)]
            # latter.clusters <- tmp.position[(latter.child.index+1):length(tmp.position)]
            former.min.dist <- min(disMat[as.character(former.child), as.character(former.clusters)])
            latter.min.dist <- min(disMat[as.character(latter.child), as.character(former.clusters)])
            if (former.min.dist > latter.min.dist) {
                # switch
                tmp.position[former.child.index] <- latter.child
                tmp.position[latter.child.index] <- former.child
                message(sprintf("switch %s %s at level %s", former.child, latter.child, tmp.level))
            }
        }
        # sort
        tmp.levels <- bt2m.cellMeta[,tmp.level]
        tmp.levels <- factor(tmp.levels, levels = tmp.position)
        bt2m.cellMeta <- bt2m.cellMeta[order(tmp.levels, decreasing = F),]
        # break
    }
    bt2m.result[["cellMeta"]] <- bt2m.cellMeta
    return(bt2m.result)
}

#' rename the clusters inside the bt2m result
#'
#' @param bt2m.result A result file from RunBt2m() function

#' @return A renamed list. cellMeta contains the final bifurcation for each level
#' marker_chain contains all the significant markers for each cluster
#' bifucation contains the bifurcation details (parent, child1, child2)
#' @export
#'
RenameBt2m <- function(bt2m.result) {
  # rename three files, bt2m.cellMeta is one type, bt2m.marker.chain and bt2m.bifucation are another type
  bt2m.cellMeta <- bt2m.result[["cellMeta"]]
  bt2m.marker.chain <- bt2m.result[["marker_chain"]]
  bt2m.bifucation <- bt2m.result[["bifucation"]]
  #
  # reoder the name of cluster name and marker module name
  all.map.df <- data.frame()
  for (i in 1:ncol(bt2m.cellMeta)) {
    # find a key error, we start from L1
    # all old cluster name in each level
    tmp.cluster.name <- unique(bt2m.cellMeta[,i])
    # add level prefix
    tmp.cluster.name2 <- paste("L",i-1,"_",tmp.cluster.name,sep = "")
    # give new cluster name
    new.cluster.name <- 1:length(tmp.cluster.name2)
    tmp.map.df <- data.frame(old_cluster=tmp.cluster.name, old_cluster_prefix=tmp.cluster.name2,
                             new_cluster=new.cluster.name, stringsAsFactors = F, level=paste("L",i-1,sep = ""))
    all.map.df <- rbind(all.map.df, tmp.map.df)
    # replace old cluster name by new cluster name
    bt2m.cellMeta[,i] <- plyr::mapvalues(bt2m.cellMeta[,i],
                                            from = tmp.map.df$old_cluster,
                                            to=tmp.map.df$new_cluster)
    # break
  }
  all.map.df$new_cluster_full <- paste(all.map.df$level, all.map.df$new_cluster, sep = "_")
  #
  all.map.df.NoEnd <- subset(all.map.df, !grepl("end", old_cluster))
  # these two files don't invovle end-node
  bt2m.marker.chain$cluster <- plyr::mapvalues(bt2m.marker.chain$cluster,
                                                      from = all.map.df.NoEnd$old_cluster_prefix,
                                                      to = all.map.df.NoEnd$new_cluster_full)
  bt2m.marker.chain$parent <- plyr::mapvalues(bt2m.marker.chain$parent,
                                                    from = all.map.df.NoEnd$old_cluster_prefix,
                                                    to = all.map.df.NoEnd$new_cluster_full)
  #
  bt2m.bifucation$parent <- plyr::mapvalues(bt2m.bifucation$parent,
                                               from = all.map.df.NoEnd$old_cluster_prefix,
                                               to = all.map.df.NoEnd$new_cluster_full)
  bt2m.bifucation$child1 <- plyr::mapvalues(bt2m.bifucation$child1,
                                               from = all.map.df.NoEnd$old_cluster_prefix,
                                               to = all.map.df.NoEnd$new_cluster_full)
  bt2m.bifucation$child2 <- plyr::mapvalues(bt2m.bifucation$child2,
                                               from = all.map.df.NoEnd$old_cluster_prefix,
                                               to = all.map.df.NoEnd$new_cluster_full)
  for (i in 1:ncol(bt2m.cellMeta)) {
    # print(i)
    bt2m.cellMeta[,i] <- as.numeric(bt2m.cellMeta[,i])
    # break
  }
  # remove duplicated levels
  bt2m.cellMeta <- bt2m.cellMeta[,!duplicated(t(bt2m.cellMeta))]
  # output
  bt2m.result[["cellMeta"]] <- bt2m.cellMeta
  bt2m.result[["marker_chain"]] <- bt2m.marker.chain
  bt2m.result[["bifucation"]] <- bt2m.bifucation
  return(bt2m.result)
}

#' transform dataframe to vector
#'
#' @param df A dataframe from Seuratobj[["attr"]]

#' @return A vector with names
#' @export
#'
DataframeToVector <- function(df) {
  # transform dataframe to vector
  tmp.anno.vector <- df[,1]
  names(tmp.anno.vector) <- rownames(df)
  return(tmp.anno.vector)
}

#' filter uniquely expressed markers by expression percentage
#'
#' @param seuratObj A Seurat object
#' @param bt2m.cellMeta bt2m.cellMeta dataframe from bt2m
#' @param bt2m.marker.chain bt2m.marker.chain dataframe contains all the markers
#' @param assay Assay used for prediction
#' @param slot slot used for prediction
#'
#' @return A new dataframe with two additional columns: cluster_pct (expression percentage in the cluster) and bcg_pct (expression percentage in the background cells)
#' @export
#'
AddMarkerExpressionPct <- function(seuratObj, bt2m.cellMeta, bt2m.marker.chain, assay = "RNA", slot = "data") {
  # # cannot remove any markers, they all meaningful
  # rm duplicate markers from last
  # bt2m.marker.chain.uniq <- RemoveDuplicatedMarker(bt2m.marker.chain)
  bt2m.marker.chain$cluster_pct <- 0
  bt2m.marker.chain$bcg_pct <- 0
  # get data
  reference.data <- GetAssayData(object = seuratObj, assay = assay, slot = slot)
  # tmp.exprMat <- as.matrix(seuratObj@assays$RNA@data)
  # too slow!!!
  # for (i in 1:nrow(bt2m.marker.chain)) {
  #   # print(i)
  #   tmp.gene <- bt2m.marker.chain$gene[i]
  #   tmp.cluster <- bt2m.marker.chain$cluster[i]
  #   tmp.level <- strsplit(tmp.cluster, split = "_")[[1]][1]
  #   tmp.cluster.index <- as.numeric(strsplit(tmp.cluster, split = "_")[[1]][2])
  #   tmp.cells <- rownames(bt2m.cellMeta[bt2m.cellMeta[,tmp.level]==tmp.cluster.index,] )
  #   bcg.cells <- rownames(bt2m.cellMeta[bt2m.cellMeta[,tmp.level]!=tmp.cluster.index,] )
  #   tmp.gene.expression <- reference.data[tmp.gene,tmp.cells]
  #   bcg.gene.expression <- reference.data[tmp.gene,bcg.cells]
  #   tmp.cluster.pct <- sum(tmp.gene.expression>0.1) / length(tmp.gene.expression)
  #   tmp.background.pct <- sum(bcg.gene.expression>0.1) / length(bcg.gene.expression)
  #   bt2m.marker.chain[i,"cluster_pct"] <- tmp.cluster.pct
  #   bt2m.marker.chain[i,"bcg_pct"] <- tmp.background.pct
  #   # break
  # }
  # to reduce the computation cost, must do it cluster by cluster
  for (tmp.cluster in unique(bt2m.marker.chain$cluster)) {
      # print(i)
      tmp.df <- subset(bt2m.marker.chain, cluster==tmp.cluster)
      tmp.level <- strsplit(tmp.cluster, split = "_")[[1]][1]
      tmp.cluster.index <- as.numeric(strsplit(tmp.cluster, split = "_")[[1]][2])
      tmp.cells <- rownames(bt2m.cellMeta[bt2m.cellMeta[,tmp.level]==tmp.cluster.index,] )
      bcg.cells <- rownames(bt2m.cellMeta[bt2m.cellMeta[,tmp.level]!=tmp.cluster.index,] )
      tmp.expression <- reference.data[,tmp.cells]
      bcg.expression <- reference.data[,bcg.cells]
      tmp.cluster.pct <- rowSums(as.matrix(tmp.expression>0)) / ncol(tmp.expression)
      tmp.background.pct <- rowSums(as.matrix(bcg.expression>0)) / ncol(bcg.expression)
      # write to df
      bt2m.marker.chain[bt2m.marker.chain$gene %in% tmp.df$gene & bt2m.marker.chain$cluster==tmp.cluster,]$cluster_pct <- tmp.cluster.pct[tmp.df$gene]
      bt2m.marker.chain[bt2m.marker.chain$gene %in% tmp.df$gene & bt2m.marker.chain$cluster==tmp.cluster,]$bcg_pct <- tmp.background.pct[tmp.df$gene]
      # break
  }
  #
  return(bt2m.marker.chain)
}

#' remove duplicated GO terms based on overlaps of genes
#'
#' @param tmp.GO.df GO annotation dataframe from clusterProfiler package
#' @param max.overlap Maximum overlapped percentage of genes
#'
#' @return GO annotation dataframe without duplications
#' @export
#'
RemoveDuplicatedGO <- function(tmp.GO.df, max.overlap=0.6) {
  # remove duplicated GO
  remain.GOterms <- c()
  # remove duplicate GO or KEGG terms
  for (i in 1:length(tmp.GO.df$ID)) {
    if (i == 1) {remain.GOterms <- c(remain.GOterms, i); next}
    leftGenes <- strsplit(paste(tmp.GO.df[remain.GOterms,]$geneID, collapse = "/"), split = "/")[[1]]
    tmpGenes <- strsplit(paste(tmp.GO.df[i,]$geneID, collapse = "/"), split = "/")[[1]]
    # if (length(tmpGenes)<2) {next}
    overlap_genes <- intersect(leftGenes, tmpGenes)
    if (length(overlap_genes)/length(tmpGenes) > 0.6) {next}
    remain.GOterms <- c(remain.GOterms, i)
  }
  return(tmp.GO.df[remain.GOterms,])
}

#' Perform GO annotation of bt2m.marker.chain dataframe
#'
#' @param bt2m.marker.chain terbi.marker.chain dataframe contains all the markers
#' @param organism "hs" for Homo sapiens, or "mm" for Mus musculus
#' @param pvalueCutoff cut off P-value for GO annotations
#' @param min_count Minimal count of genes for GO annotations
#'
#' @return GO annotation dataframe labeled with cluster
#' @export
#'
Bt2mEnrichGO <- function(bt2m.marker.chain, organism="hs", pvalueCutoff = 0.05, min_count = 3) {
  marker.list <- list()
  for (i in unique(bt2m.marker.chain$cluster)) {
    # print(i)
    tmp.marker.df <- subset(bt2m.marker.chain, cluster==i)
    marker.list[[i]] <- tmp.marker.df$gene
  }
  #
  bt2m.GO.anno <- clusterProfilerORA(geneList = marker.list, organism=organism)
  #
  bt2m.GO.anno.filtered <- data.frame()
  for (i in names(bt2m.GO.anno$go_list)) {
    # print(i)
    tmp.GO.df <- bt2m.GO.anno$go_list[[i]]@result
    tmp.GO.df <- subset(tmp.GO.df, pvalue<pvalueCutoff & p.adjust<pvalueCutoff & qvalue<pvalueCutoff & Count>min_count)
    if (nrow(tmp.GO.df)<1) next
    tmp.GO.df <- RemoveDuplicatedGO(tmp.GO.df)
    tmp.GO.df$cluster <- i
    bt2m.GO.anno.filtered <- rbind(bt2m.GO.anno.filtered, tmp.GO.df)
    # break
  }
  return(bt2m.GO.anno.filtered)
}

#' GO KEGG ORA analysis by clusterProfiler
#' @param geneList a gene list
#' @param organism "hs" for Homo sapiens, or "mm" for Mus musculus
#'
#' @return a list with GO and KEGG annotation result (clusterProfiler format)
#' @export
#'
clusterProfilerORA <- function(geneList=markerList, organism="hs") {
  # library(clusterProfiler)
  go_list <- list()
  kegg_list <- list()
  nameList <- names(geneList)
  if (is.null(nameList)) {
    print("no name for the gene list!!!\n")
    nameList <- 1:length(geneList)
  }
  for (i in nameList) {
    genes <- geneList[[i]]
    projectName <- i
    if (organism=="mm") {
      # library(org.Mm.eg.db) # mouse
      gene.df <- bitr(genes, fromType = "SYMBOL", toType = c("ENSEMBL", "ENTREZID"), OrgDb = org.Mm.eg.db)
      ego <- enrichGO(gene      = gene.df$ENTREZID,
                  #universe      = genes, # SYMBOL
                  keyType       = "ENTREZID",
                  OrgDb         = org.Mm.eg.db,
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05,
                  readable      = TRUE)
      # remove duplications # too slow
      # ego <- clusterProfiler::simplify(ego, cutoff=0.7, by="p.adjust", select_fun=min)
      # drop top 3 level
      # ego <- dropGO(ego, level = c(1,2))
      go_list[[i]] <- ego
      kk <- enrichKEGG(gene = gene.df$ENTREZID, organism = 'mmu', pvalueCutoff = 0.05)
      # kk$genes <- unlist(lapply(kk$geneID, ID2gene))
      kegg_list[[i]] <- kk
      }
    else if (organism=="hs") {
      # library(org.Hs.eg.db)
      gene.df <- bitr(genes, fromType = "SYMBOL", toType = c("ENSEMBL", "ENTREZID"), OrgDb = org.Hs.eg.db)
      ego <- enrichGO(gene      = gene.df$ENTREZID,
                  #universe      = genes, # SYMBOL
                  keyType       = "ENTREZID",
                  OrgDb         = org.Hs.eg.db,
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05,
                  readable      = TRUE)
      # remove duplications
      # too slow
      # ego <- clusterProfiler::simplify(ego, cutoff=0.7, by="p.adjust", select_fun=min)
      # drop top 3 level
      # ego <- dropGO(ego, level = c(1,2))
      go_list[[i]] <- ego
      kk <- enrichKEGG(gene = gene.df$ENTREZID, organism = 'hsa', pvalueCutoff = 0.05)
      # kk$genes <- unlist(lapply(kk$geneID, ID2gene))
      kegg_list[[i]] <- kk
    }
    else {stop("only support hs and mm now!")}
  }
    return(list("go_list"=go_list, "kegg_list"=kegg_list))
}

#' Subset seurat object by cell names
#'
#' @param seuratObj A Seurat object
#' @param tmp.cells cells for subseting
#'
#' @return A subset of seurat object
#' @export
#'
subsetSeuratObjByCells <- function(seuratObj, tmp.cells) {
  # tmp.cells <- rownames(tmp.bt2m.cellMeta.3)
  tmp <- rep(F, ncol(seuratObj))
  names(tmp) <- colnames(seuratObj)
  tmp[tmp.cells] <- T
  seuratObj$select_cells <- factor(tmp, levels = c(T, F))
  tmp.seuratObj <- subset(seuratObj, select_cells == TRUE)
  if (sum(colnames(tmp.seuratObj)!=names((tmp.seuratObj@active.ident))) != 0) stop("colnames in seuratObj and active.ident not match")
  return(tmp.seuratObj)
}

#' Find initial bifurcation event of any two clusters in bt2m.cellMeta
#'
#' @param seuratObj A Seurat object
#' @param bt2m.cellMeta bt2m.cellMeta from bt2m
#' @param cluster_1 cluster 1
#' @param cluster_2 cluster 2
#' @param balance_cells output balanced cell number (T or F)
#'
#' @return A list. subset_seuratObj is the subset cells in cluster 1 and cluster 2
#' diff_marker_1 is the marker of cluster 1, diff_marker_2 is the marker of cluster 2
#' @export
#'
GetInitialBifurcation <- function(seuratObj, bt2m.cellMeta, cluster_1, cluster_2, balance_cells=T) {
  # input any two cluster, trace back to the initial bifurcation event
  level_1 <- strsplit(cluster_1, split = "_")[[1]][1]
  cluster_1.index <- as.numeric(strsplit(cluster_1, split = "_")[[1]][2])
  level_1.index <- as.numeric(substr(level_1, 2, 10))
  #
  level_2 <- strsplit(cluster_2, split = "_")[[1]][1]
  cluster_2.index <- as.numeric(strsplit(cluster_2, split = "_")[[1]][2])
  level_2.index <- as.numeric(substr(level_2, 2, 10))
  # seperate bt2m.cellMeta file
  tmp.bt2m.cellMeta.1 <- bt2m.cellMeta[bt2m.cellMeta[,level_1]==cluster_1.index,]
  tmp.bt2m.cellMeta.2 <- bt2m.cellMeta[bt2m.cellMeta[,level_2]==cluster_2.index,]
  #
  # if they are inherited relationship
  if (!(row.names(tmp.bt2m.cellMeta.1) %in% row.names(tmp.bt2m.cellMeta.2))==0 || !(row.names(tmp.bt2m.cellMeta.2) %in% row.names(tmp.bt2m.cellMeta.1))==0) {
    stop("Input clusters are inherited. Please select two unrelated clusters!")
  } else {
    # 倒序搜索共同祖先
    #tmp.bt2m.cellMeta.3 <- bt2m.cellMeta[bt2m.cellMeta[,level_1]==cluster_1.index | bt2m.cellMeta[,level_2]==cluster_2.index,]
    tmp.bt2m.cellMeta.3 <- rbind(tmp.bt2m.cellMeta.1, tmp.bt2m.cellMeta.2)
    for (i in min(level_1.index, level_2.index):1) {
      tmp.level <- paste("L", i-1, sep = "")
      # print(tmp.level)
      if (length(unique(tmp.bt2m.cellMeta.3[,tmp.level]))==1) break
    }
    None.bifucation.level.index <- i
    None.bifucation.level <- paste("L", None.bifucation.level.index-1, sep = "")
    alreay.bifucation.level <- paste("L", None.bifucation.level.index, sep = "")
    # get the first two ancestor
    ancestor.cluster_1 <- tmp.bt2m.cellMeta.3[tmp.bt2m.cellMeta.3[,level_1]==cluster_1.index,][1,alreay.bifucation.level]
    ancestor.cluster_2 <- tmp.bt2m.cellMeta.3[tmp.bt2m.cellMeta.3[,level_2]==cluster_2.index,][1,alreay.bifucation.level]
    #
    Non.split.ancestor.cluster <- tmp.bt2m.cellMeta.3[tmp.bt2m.cellMeta.3[,level_1]==cluster_1.index,][1,None.bifucation.level]
    # get the first bifurcation markers
    diff_marker_1 <- subset(bt2m.marker.chain, cluster==paste(alreay.bifucation.level, ancestor.cluster_1, sep = "_"))
    diff_marker_2 <- subset(bt2m.marker.chain, cluster==paste(alreay.bifucation.level, ancestor.cluster_2, sep = "_"))
    # get subset seuratObj
    tmp.seuratObj <- subsetSeuratObj(seuratObj, rownames(tmp.bt2m.cellMeta.3))
    tmp.bt2m.cellMeta.3$merged_cluster <- paste(level_2, tmp.bt2m.cellMeta.3[,level_2], sep = "_")
    tmp.bt2m.cellMeta.3[tmp.bt2m.cellMeta.3[,level_1]==cluster_1.index, "merged_cluster"] <- paste(level_1, tmp.bt2m.cellMeta.3[tmp.bt2m.cellMeta.3[,level_1]==cluster_1.index, level_1], sep = "_")
    # set merged cluster name
    tmp.ident <- tmp.bt2m.cellMeta.3[colnames(tmp.seuratObj),]$merged_cluster
    names(tmp.ident) <- colnames(tmp.seuratObj)
    tmp.ident <- factor(tmp.ident, levels = c(cluster_1, cluster_2))
    tmp.seuratObj@active.ident <- tmp.ident
    #
    if (balance_cells) {
      # balance cell number
      set.seed(49)
      cell_number_table <- table(tmp.seuratObj@active.ident)
      cell.1 <- sample(colnames(tmp.seuratObj)[as.character(tmp.seuratObj@active.ident)==names(cell_number_table)[1]], size = min(cell_number_table))
      cell.2 <- sample(colnames(tmp.seuratObj)[as.character(tmp.seuratObj@active.ident)==names(cell_number_table)[2]], size = min(cell_number_table))
      tmp.seuratObj <- subsetSeuratObj(tmp.seuratObj, c(cell.1, cell.2))

    }
    return(list(subset_seuratObj=tmp.seuratObj, diff_marker_1=diff_marker_1, diff_marker_2=diff_marker_2))
  }
}

#' Find the cluster chain of any given cluster
#'
#' @param bt2m.cellMeta bt2m.cellMeta from bt2m
#' @param bt2m.bifucation bt2m.bifucation from bt2m
#' @param target.cluster target cluster
#'
#' @return A vector, names are cluster name, element shows the states of clusters
#' @export
#'
GetClusterChain <- function(bt2m.cellMeta, bt2m.bifucation, target.cluster) {
    tmp.level <- strsplit(target.cluster, split = "_")[[1]][1]
    tmp.level.index <- as.numeric(substr(tmp.level, 2, 10))
    tmp.cluster.index <- as.numeric(strsplit(target.cluster, split = "_")[[1]][2])
    if (!(tmp.level %in% colnames(bt2m.cellMeta)) || !(tmp.cluster.index %in% bt2m.cellMeta[,tmp.level])) {
        message("Please input correct target.cluster!!!")
    } else {
        cluster.chain <- paste(colnames(bt2m.cellMeta[bt2m.cellMeta[,tmp.level]==tmp.cluster.index,]), 
                               bt2m.cellMeta[bt2m.cellMeta[,tmp.level]==tmp.cluster.index,][1,], sep = "_")
    }
    # node attribute
    node_attr <- c("End_node","Parent_node")[as.integer(cluster.chain %in% unlist(bt2m.bifucation))+1]
    names(node_attr) <- cluster.chain
    return(node_attr)
}

#' Write bt2m result to Seurat object
#'
#' @param seuratObj A Seurat object
#' @param bt2m.result bt2m.result
#'
#' @return A Seurat object contains bt2m.result. bt2m will be stored in assay data region of Seurat object ("seuratObj@assays$bt2m")
#' @export
#'
WriteBt2mIntoSeurat <- function(seuratObj, bt2m.result) {
  # get results
  bt2m.cellMeta <- bt2m.result[["cellMeta"]]
  bt2m.marker.chain <- bt2m.result[["marker_chain"]]
  bt2m.bifucation <- bt2m.result[["bifucation"]]
  # input L-data to seurat
  for (i in colnames(bt2m.cellMeta)) {
    # print(i)
    tmp.ident <- bt2m.cellMeta[,i]
    names(tmp.ident) <- rownames(bt2m.cellMeta)
    tmp.ident <- factor(tmp.ident, levels = sort(unique(tmp.ident)))
    seuratObj[[i]] <- tmp.ident
  }
  #bt2m.cellMeta <- bt2m.result[["cellMeta"]]
  #bt2m.marker.chain <- bt2m.result[["marker_chain"]]
  #bt2m.bifucation <- bt2m.result[["bifucation"]]
  seuratObj@assays$bt2m <- bt2m.result
  return(seuratObj)
}

#' Remove duplicated markers (just for plotting)
#'
#' @param bt2m.marker.chain bt2m.marker.chain from bt2m
#' @param method select a method to sort markers
#'
#' @return bt2m.marker.chain without duplicated genes
#' @export
#'
RemoveDuplicatedMarker <- function(bt2m.marker.chain, method = "correlation") {
  if (method == "level") {
    bt2m.marker.chain.uniq <- bt2m.marker.chain[!duplicated(bt2m.marker.chain$gene, fromLast = T),]
    } else if (method == "correlation") {
      bt2m.marker.chain <- bt2m.marker.chain[order(bt2m.marker.chain$correlation, decreasing = T),]
      bt2m.marker.chain.uniq <- bt2m.marker.chain[!duplicated(bt2m.marker.chain$gene, fromLast = F),]
    } else {
      message("Please input correct method!")
    }
  return(bt2m.marker.chain.uniq)
}
