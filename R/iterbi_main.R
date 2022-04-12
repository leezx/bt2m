#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
globalVariables(
  names=c("Count","Description","avg_logFC","bcg_pct","cluster","cluster_pct","correlation","iterbi.cellMeta","iterbi.marker.chain","markerList","new_cluster","old_cluster","p_val","p_val_adj","pvalue","qvalue","select_cells","seuratObj"),
  package = "iterbi",
  add = TRUE
)
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
RunIterbi <- function(seuratObj, method = "graph", min.marker.num = 100, max.level.num = 20,
                       min.cell.count = 50, resolution.sets = 30, verbose = T) {
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
  iterbi.cellMeta <- data.frame(row.names = seuratObj$cellName, cellName=seuratObj$cellName, stringsAsFactors = F)
  # add initial one
  iterbi.cellMeta[,"L0"] <- 1
  # add default L1-max.level.num
  for (i in paste("L", 1:max.level.num, sep = "")) {
    # print(i)
    iterbi.cellMeta[,i] <- 0
  }
  #
  # cannot start from specifc level, the marker list will be chaos
  iterbi.marker.chain <- data.frame()
  iterbi.bifucation <- data.frame()
  # iteratively bifurcation until no enough markers can be found
  for (i in 0:max.level.num) {
    # prevoius.level.index <- paste("L",i-1,sep = "")
    tmp.level.index <- paste("L",i,sep = "") # for getting data
    next.level.index <- paste("L",i+1,sep = "") # for writing data
    message(paste("We are now at", tmp.level.index, sep = " "))
    # small testing
    # if (tmp.level.index=="L3") break
    #
    for (j in unique(iterbi.cellMeta[,tmp.level.index])) {
      # print(j)
      tmp.cluster.index <- paste(tmp.level.index, j, sep = "_")
      # get cells
      # message("get tmp.cells...")
      tmp.cells <- iterbi.cellMeta[iterbi.cellMeta[,tmp.level.index]==j,]$cellName
      # if there are no marker for the previous level, endpoint
      if (grepl("end", j)) {
        iterbi.cellMeta[tmp.cells, next.level.index] <- j
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
        iterbi.cellMeta[tmp.cells, next.level.index] <- paste("end", next.level.index, 2*j-1, sep = "-")
        # message(paste(tmp.cluster.index, "has no enough cells, this is a end point...", sep = " "))
        next
      }
      #
      # binary clustering
      # message(length(tmp.cells))
      # subset funciton have problem
      tmp <- rep(F, ncol(seuratObj))
      names(tmp) <- colnames(seuratObj)
      tmp[tmp.cells] <- T
      seuratObj$select_cells <- factor(tmp, levels = c(T, F))
      tmp.seuratObj <- subset(seuratObj, subset = select_cells == T)
      # tmp.seuratObj <- subset(seuratObj, subset = cellName %in% tmp.cells) # can't find tmp.cells, don't know why?
      #
      if (method == "graph") {tmp.seuratObj <- IterbiBifucation.graph(tmp.seuratObj, resolution.sets = resolution.sets)}
      else if (method == "hclust") {tmp.seuratObj <- IterbiBifucation.hclust(tmp.seuratObj)}
      else if (method == "kmeans") {tmp.seuratObj <- IterbiBifucation.kmeans(tmp.seuratObj)}
      else {message("Please select one method in: graph, hclust, kmeans!")}
      # check clustering result
      if (length(unique(tmp.seuratObj@active.ident)) != 2) {
        message("Cannot do bifurcation. Set bigger resolution.sets!!!")
        stop()
      }
      # assign name to the marker dfs
      next.cluster.index.1 <- paste(next.level.index, 2*j-1, sep = "_")
      next.cluster.index.2 <- paste(next.level.index, 2*j, sep = "_")
      # write annotation to the next level
      iterbi.cellMeta[tmp.seuratObj$cellName, next.level.index] <- as.integer(tmp.seuratObj@active.ident) + 2*(j-1)
      # get markers
      binary_markers <- IdentifyBinaryMarkers(tmp.seuratObj)
      b1.markers <- binary_markers[["b1"]]
      b2.markers <- binary_markers[["b2"]]
      # set end point
      if (nrow(b1.markers) >= min.marker.num) {
        b1.markers$parent <- tmp.cluster.index
        b1.markers$cluster <- next.cluster.index.1
        iterbi.marker.chain <- rbind(iterbi.marker.chain, b1.markers)
      }
      if (nrow(b2.markers) >= min.marker.num) {
        b2.markers$parent <- tmp.cluster.index
        b2.markers$cluster <- next.cluster.index.2
        iterbi.marker.chain <- rbind(iterbi.marker.chain, b2.markers)
      }
      # final judgement
      if ((nrow(b1.markers) < min.marker.num) && (nrow(b2.markers) < min.marker.num)) {
        # message(paste(tmp.cluster.index, "has no marker, this is a end point...", sep = " "))
        # set end point (take care, what the name is, when no marker)
        # must add "end" and "level" to avoid the mixture between levels
        iterbi.cellMeta[tmp.cells, next.level.index] <- paste("end", next.level.index, 2*j-1, sep = "-")
      } else {
        tmp.bifucation <- data.frame(parent=tmp.cluster.index, child1=next.cluster.index.1, child2=next.cluster.index.2)
        iterbi.bifucation <- rbind(iterbi.bifucation, tmp.bifucation)
        if (verbose) message(sprintf("Successfully split %s to %s and %s", tmp.cluster.index,
                                     next.cluster.index.1, next.cluster.index.2))
      }
    }
    # exist running when all cluster are end points
    if (sum(!grepl("end", unique(iterbi.cellMeta[,tmp.level.index]))) == 0) {
      message("Bifurcating stopped! No more clusters can be split")
      break
    }
  }
  # remove emplty levels
  iterbi.cellMeta <- iterbi.cellMeta[,colSums(iterbi.cellMeta!=0)!=0]
  iterbi.cellMeta$cellName <- NULL
  # remove same columns
  iterbi.cellMeta <- iterbi.cellMeta[,!duplicated(t(iterbi.cellMeta))]
  # sort the df last to first
  for (i in ncol(iterbi.cellMeta):1) {
    iterbi.cellMeta <- iterbi.cellMeta[order( iterbi.cellMeta[,i] ),]
  }
  return(list(cellMeta=iterbi.cellMeta, marker_chain=iterbi.marker.chain, bifucation=iterbi.bifucation))
}

#' rename the clusters inside the iterbi result
#'
#' @param iterbi.result A result file from RunIterbi() function

#' @return A renamed list. cellMeta contains the final bifurcation for each level
#' marker_chain contains all the significant markers for each cluster
#' bifucation contains the bifurcation details (parent, child1, child2)
#' @export
#'
RenameIterbi <- function(iterbi.result) {
  # rename three files, iterbi.cellMeta is one type, iterbi.marker.chain and iterbi.bifucation are another type
  iterbi.cellMeta <- iterbi.result[["cellMeta"]]
  iterbi.marker.chain <- iterbi.result[["marker_chain"]]
  iterbi.bifucation <- iterbi.result[["bifucation"]]
  #
  # reoder the name of cluster name and marker module name
  all.map.df <- data.frame()
  for (i in 1:ncol(iterbi.cellMeta)) {
    # find a key error, we start from L1
    # all old cluster name in each level
    tmp.cluster.name <- unique(iterbi.cellMeta[,i])
    # add level prefix
    tmp.cluster.name2 <- paste("L",i-1,"_",tmp.cluster.name,sep = "")
    # give new cluster name
    new.cluster.name <- 1:length(tmp.cluster.name2)
    tmp.map.df <- data.frame(old_cluster=tmp.cluster.name, old_cluster_prefix=tmp.cluster.name2,
                             new_cluster=new.cluster.name, stringsAsFactors = F, level=paste("L",i-1,sep = ""))
    all.map.df <- rbind(all.map.df, tmp.map.df)
    # replace old cluster name by new cluster name
    iterbi.cellMeta[,i] <- plyr::mapvalues(iterbi.cellMeta[,i],
                                            from = tmp.map.df$old_cluster,
                                            to=tmp.map.df$new_cluster)
    # break
  }
  all.map.df$new_cluster_full <- paste(all.map.df$level, all.map.df$new_cluster, sep = "_")
  #
  all.map.df.NoEnd <- subset(all.map.df, !grepl("end", old_cluster))
  # these two files don't invovle end-node
  iterbi.marker.chain$cluster <- plyr::mapvalues(iterbi.marker.chain$cluster,
                                                      from = all.map.df.NoEnd$old_cluster_prefix,
                                                      to = all.map.df.NoEnd$new_cluster_full)
  iterbi.marker.chain$parent <- plyr::mapvalues(iterbi.marker.chain$parent,
                                                    from = all.map.df.NoEnd$old_cluster_prefix,
                                                    to = all.map.df.NoEnd$new_cluster_full)
  #
  iterbi.bifucation$parent <- plyr::mapvalues(iterbi.bifucation$parent,
                                               from = all.map.df.NoEnd$old_cluster_prefix,
                                               to = all.map.df.NoEnd$new_cluster_full)
  iterbi.bifucation$child1 <- plyr::mapvalues(iterbi.bifucation$child1,
                                               from = all.map.df.NoEnd$old_cluster_prefix,
                                               to = all.map.df.NoEnd$new_cluster_full)
  iterbi.bifucation$child2 <- plyr::mapvalues(iterbi.bifucation$child2,
                                               from = all.map.df.NoEnd$old_cluster_prefix,
                                               to = all.map.df.NoEnd$new_cluster_full)
  for (i in 1:ncol(iterbi.cellMeta)) {
    # print(i)
    iterbi.cellMeta[,i] <- as.numeric(iterbi.cellMeta[,i])
    # break
  }
  # remove duplicated levels
  iterbi.cellMeta <- iterbi.cellMeta[,!duplicated(t(iterbi.cellMeta))]
  # output
  iterbi.result[["cellMeta"]] <- iterbi.cellMeta
  iterbi.result[["marker_chain"]] <- iterbi.marker.chain
  iterbi.result[["bifucation"]] <- iterbi.bifucation
  return(iterbi.result)
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
#' @param iterbi.marker.chain iterbi.marker.chain dataframe contains all the markers
#' @param assay Assay used for prediction
#' @param slot slot used for prediction
#' @param min.cluster.pct Minimal expression percentage of target cluster
#' @param max.bcg.pct Maximum expression percentage of the background cells (non-target cells)
#' @param min.diff Minimal difference (expression percentage) between target cluster and background cells
#'
#' @return Uniquely expressed iterbi.marker.chain
#' @export
#'
GetUniqueMarker <- function(seuratObj, iterbi.marker.chain, assay = "RNA", slot = "data",
                              min.cluster.pct = 0.3, max.bcg.pct = 0.1, min.diff = 0.3) {
  # rm duplicate markers from last
  iterbi.marker.chain.uniq <- RemoveDuplicatedMarker(iterbi.marker.chain)
  iterbi.marker.chain.uniq$cluster_pct <- 0
  iterbi.marker.chain.uniq$bcg_pct <- 0
  # get data
  reference.data <- GetAssayData(object = seuratObj, assay = assay, slot = slot)
  # tmp.exprMat <- as.matrix(seuratObj@assays$RNA@data)
  for (i in 1:nrow(iterbi.marker.chain.uniq)) {
    # print(i)
    tmp.gene <- iterbi.marker.chain.uniq$gene[i]
    tmp.cluster <- iterbi.marker.chain.uniq$cluster[i]
    tmp.level <- strsplit(tmp.cluster, split = "_")[[1]][1]
    tmp.cluster.index <- as.numeric(strsplit(tmp.cluster, split = "_")[[1]][2])
    tmp.cells <- rownames(iterbi.cellMeta[iterbi.cellMeta[,tmp.level]==tmp.cluster.index,] )
    bcg.cells <- rownames(iterbi.cellMeta[iterbi.cellMeta[,tmp.level]!=tmp.cluster.index,] )
    tmp.gene.expression <- reference.data[tmp.gene,tmp.cells]
    bcg.gene.expression <- reference.data[tmp.gene,bcg.cells]
    tmp.cluster.pct <- sum(tmp.gene.expression>0.1) / length(tmp.gene.expression)
    tmp.background.pct <- sum(bcg.gene.expression>0.1) / length(bcg.gene.expression)
    iterbi.marker.chain.uniq[i,"cluster_pct"] <- tmp.cluster.pct
    iterbi.marker.chain.uniq[i,"bcg_pct"] <- tmp.background.pct
    # break
  }
  # filter
  iterbi.marker.chain.uniq <- subset(iterbi.marker.chain.uniq, (cluster_pct-bcg_pct)>0.3 && cluster_pct>min.cluster.pct && bcg_pct<max.bcg.pct)
  #
  return(iterbi.marker.chain.uniq)
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

#' Perform GO annotation of iterbi.marker.chain dataframe
#'
#' @param iterbi.marker.chain terbi.marker.chain dataframe contains all the markers
#' @param organism "hs" for Homo sapiens, or "mm" for Mus musculus
#' @param pvalueCutoff cut off P-value for GO annotations
#' @param min_count Minimal count of genes for GO annotations
#'
#' @return GO annotation dataframe labeled with cluster
#' @export
#'
IterbiEnrichGO <- function(iterbi.marker.chain, organism="hs", pvalueCutoff = 0.05, min_count = 5) {
  marker.list <- list()
  for (i in unique(iterbi.marker.chain$cluster)) {
    # print(i)
    tmp.marker.df <- subset(iterbi.marker.chain, cluster==i)
    marker.list[[i]] <- tmp.marker.df$gene
  }
  #
  iterbi.GO.anno <- clusterProfilerORA(geneList = marker.list, organism=organism)
  #
  iterbi.GO.anno.filtered <- data.frame()
  for (i in names(iterbi.GO.anno$go_list)) {
    # print(i)
    tmp.GO.df <- iterbi.GO.anno$go_list[[i]]@result
    tmp.GO.df <- subset(tmp.GO.df, pvalue<pvalueCutoff & p.adjust<pvalueCutoff & qvalue<pvalueCutoff & Count>min_count)
    if (nrow(tmp.GO.df)<1) next
    tmp.GO.df <- RemoveDuplicatedGO(tmp.GO.df)
    tmp.GO.df$cluster <- i
    iterbi.GO.anno.filtered <- rbind(iterbi.GO.anno.filtered, tmp.GO.df)
    # break
  }
  return(iterbi.GO.anno.filtered)
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
subsetSeuratObj <- function(seuratObj, tmp.cells) {
  # tmp.cells <- rownames(tmp.iterbi.cellMeta.3)
  tmp <- rep(F, ncol(seuratObj))
  names(tmp) <- colnames(seuratObj)
  tmp[tmp.cells] <- T
  seuratObj$select_cells <- factor(tmp, levels = c(T, F))
  tmp.seuratObj <- subset(seuratObj, subset = select_cells == T)
  if (sum(colnames(tmp.seuratObj)!=names((tmp.seuratObj@active.ident))) != 0) stop("colnames in seuratObj and active.ident not match")
  return(tmp.seuratObj)
}

#' Find initial bifurcation event of any two clusters in iterbi.cellMeta
#'
#' @param seuratObj A Seurat object
#' @param iterbi.cellMeta iterbi.cellMeta from iterbi
#' @param cluster_1 cluster 1
#' @param cluster_2 cluster 2
#' @param balance_cells output balanced cell number (T or F)
#'
#' @return A list. subset_seuratObj is the subset cells in cluster 1 and cluster 2
#' diff_marker_1 is the marker of cluster 1, diff_marker_2 is the marker of cluster 2
#' @export
#'
GetInitialBifurcation <- function(seuratObj, iterbi.cellMeta, cluster_1, cluster_2, balance_cells=T) {
  # input any two cluster, trace back to the initial bifurcation event
  level_1 <- strsplit(cluster_1, split = "_")[[1]][1]
  cluster_1.index <- as.numeric(strsplit(cluster_1, split = "_")[[1]][2])
  level_1.index <- as.numeric(substr(level_1, 2, 10))
  #
  level_2 <- strsplit(cluster_2, split = "_")[[1]][1]
  cluster_2.index <- as.numeric(strsplit(cluster_2, split = "_")[[1]][2])
  level_2.index <- as.numeric(substr(level_2, 2, 10))
  # seperate iterbi.cellMeta file
  tmp.iterbi.cellMeta.1 <- iterbi.cellMeta[iterbi.cellMeta[,level_1]==cluster_1.index,]
  tmp.iterbi.cellMeta.2 <- iterbi.cellMeta[iterbi.cellMeta[,level_2]==cluster_2.index,]
  #
  # if they are inherited relationship
  if (!(row.names(tmp.iterbi.cellMeta.1) %in% row.names(tmp.iterbi.cellMeta.2))==0 || !(row.names(tmp.iterbi.cellMeta.2) %in% row.names(tmp.iterbi.cellMeta.1))==0) {
    stop("Input clusters are inherited. Please select two unrelated clusters!")
  } else {
    # 倒序搜索共同祖先
    #tmp.iterbi.cellMeta.3 <- iterbi.cellMeta[iterbi.cellMeta[,level_1]==cluster_1.index | iterbi.cellMeta[,level_2]==cluster_2.index,]
    tmp.iterbi.cellMeta.3 <- rbind(tmp.iterbi.cellMeta.1, tmp.iterbi.cellMeta.2)
    for (i in min(level_1.index, level_2.index):1) {
      tmp.level <- paste("L", i-1, sep = "")
      # print(tmp.level)
      if (length(unique(tmp.iterbi.cellMeta.3[,tmp.level]))==1) break
    }
    None.bifucation.level.index <- i
    None.bifucation.level <- paste("L", None.bifucation.level.index-1, sep = "")
    alreay.bifucation.level <- paste("L", None.bifucation.level.index, sep = "")
    # get the first two ancestor
    ancestor.cluster_1 <- tmp.iterbi.cellMeta.3[tmp.iterbi.cellMeta.3[,level_1]==cluster_1.index,][1,alreay.bifucation.level]
    ancestor.cluster_2 <- tmp.iterbi.cellMeta.3[tmp.iterbi.cellMeta.3[,level_2]==cluster_2.index,][1,alreay.bifucation.level]
    #
    Non.split.ancestor.cluster <- tmp.iterbi.cellMeta.3[tmp.iterbi.cellMeta.3[,level_1]==cluster_1.index,][1,None.bifucation.level]
    # get the first bifurcation markers
    diff_marker_1 <- subset(iterbi.marker.chain, cluster==paste(alreay.bifucation.level, ancestor.cluster_1, sep = "_"))
    diff_marker_2 <- subset(iterbi.marker.chain, cluster==paste(alreay.bifucation.level, ancestor.cluster_2, sep = "_"))
    # get subset seuratObj
    tmp.seuratObj <- subsetSeuratObj(seuratObj, rownames(tmp.iterbi.cellMeta.3))
    tmp.iterbi.cellMeta.3$merged_cluster <- paste(level_2, tmp.iterbi.cellMeta.3[,level_2], sep = "_")
    tmp.iterbi.cellMeta.3[tmp.iterbi.cellMeta.3[,level_1]==cluster_1.index, "merged_cluster"] <- paste(level_1, tmp.iterbi.cellMeta.3[tmp.iterbi.cellMeta.3[,level_1]==cluster_1.index, level_1], sep = "_")
    # set merged cluster name
    tmp.ident <- tmp.iterbi.cellMeta.3[colnames(tmp.seuratObj),]$merged_cluster
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
#' @param iterbi.cellMeta iterbi.cellMeta from iterbi
#' @param iterbi.bifucation iterbi.bifucation from iterbi
#' @param target.cluster target cluster
#'
#' @return A vector, names are cluster name, element shows the states of clusters
#' @export
#'
GetClusterChain <- function(iterbi.cellMeta, iterbi.bifucation, target.cluster) {
    tmp.level <- strsplit(target.cluster, split = "_")[[1]][1]
    tmp.level.index <- as.numeric(substr(tmp.level, 2, 10))
    tmp.cluster.index <- as.numeric(strsplit(target.cluster, split = "_")[[1]][2])
    if (!(tmp.level %in% colnames(iterbi.cellMeta)) || !(tmp.cluster.index %in% iterbi.cellMeta[,tmp.level])) {
        message("Please input correct target.cluster!!!")
    } else {
        cluster.chain <- paste(colnames(iterbi.cellMeta[iterbi.cellMeta[,tmp.level]==tmp.cluster.index,]), 
                               iterbi.cellMeta[iterbi.cellMeta[,tmp.level]==tmp.cluster.index,][1,], sep = "_")
    }
    # node attribute
    node_attr <- c("End_node","Parent_node")[as.integer(cluster.chain %in% unlist(iterbi.bifucation))+1]
    names(node_attr) <- cluster.chain
    return(node_attr)
}

#' Write iterbi result to Seurat object
#'
#' @param seuratObj A Seurat object
#' @param iterbi.result iterbi.result
#'
#' @return A Seurat object contains iterbi.result. iterbi will be stored in assay data region of Seurat object ("seuratObj@assays$iterbi")
#' @export
#'
WriteIterbiIntoSeurat <- function(seuratObj, iterbi.result) {
  # get results
  iterbi.cellMeta <- iterbi.result[["cellMeta"]]
  iterbi.marker.chain <- iterbi.result[["marker_chain"]]
  iterbi.bifucation <- iterbi.result[["bifucation"]]
  # input L-data to seurat
  for (i in colnames(iterbi.cellMeta)) {
    # print(i)
    tmp.ident <- iterbi.cellMeta[,i]
    names(tmp.ident) <- rownames(iterbi.cellMeta)
    tmp.ident <- factor(tmp.ident, levels = sort(unique(tmp.ident)))
    seuratObj[[i]] <- tmp.ident
  }
  #iterbi.cellMeta <- iterbi.result[["cellMeta"]]
  #iterbi.marker.chain <- iterbi.result[["marker_chain"]]
  #iterbi.bifucation <- iterbi.result[["bifucation"]]
  seuratObj@assays$iterbi <- iterbi.result
  return(seuratObj)
}

#' Remove duplicated markers from last level
#'
#' @param iterbi.marker.chain iterbi.marker.chain from iterbi
#'
#' @return iterbi.marker.chain without duplicated genes
#' @export
#'
RemoveDuplicatedMarker <- function(iterbi.marker.chain, method = "correlation") {
  if (method == "level") {
    iterbi.marker.chain.uniq <- iterbi.marker.chain[!duplicated(iterbi.marker.chain$gene, fromLast = T),]
    } else if (method == "correlation") {
      iterbi.marker.chain <- iterbi.marker.chain[order(iterbi.marker.chain$correlation, decreasing = T),]
      iterbi.marker.chain.uniq <- iterbi.marker.chain[!duplicated(iterbi.marker.chain$gene, fromLast = F),]
    } else {
      message("Please input correct method!")
    }
  return(iterbi.marker.chain.uniq)
}
