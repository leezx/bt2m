run.iterbi <- function(seuset, method = "graph", min.marker.num = 100, max.level.num = 20,
                       min.cell.count = 50, res_sets = 30, verbose = T) {
  # method = "graph"; min.marker.num = 100; max.level.num = 20; min.cell.count = 50; verbose = T # for quite test
  # suppress warnings
  options(warn = -1, stringsAsFactors = F)
  # min.marker.num <- 100
  # max.level.num <- 20
  # min.cell.count <- 50 # must be more than number PCs
  # verbose <- T
  # key file 1: the cell annotation dataframe
  seuset$cellName <- as.character(colnames(seuset))
  iterbi.cellMeta <- data.frame(row.names = seuset$cellName, cellName=seuset$cellName, stringsAsFactors = F)
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
      tmp <- rep(F, ncol(seuset))
      names(tmp) <- colnames(seuset)
      tmp[tmp.cells] <- T
      seuset$select_cells <- factor(tmp, levels = c(T, F))
      tmp.seuset <- subset(seuset, subset = select_cells == T)
      # tmp.seuset <- subset(seuset, subset = cellName %in% tmp.cells) # can't find tmp.cells, don't know why?
      #
      if (method == "graph") {tmp.seuset <- iterbi.bifucation.graph(tmp.seuset, res_sets = res_sets)}
      else if (method == "hclust") {tmp.seuset <- iterbi.bifucation.hclust(tmp.seuset)}
      else if (method == "kmeans") {tmp.seuset <- iterbi.bifucation.kmeans(tmp.seuset)}
      else {message("Please select one method in: graph, hclust, kmeans!")}
      # check clustering result
      if (length(unique(tmp.seuset@active.ident)) != 2) {
        message("Cannot do bifurcation. Set bigger res_sets!!!")
        stop()
      }
      # assign name to the marker dfs
      next.cluster.index.1 <- paste(next.level.index, 2*j-1, sep = "_")
      next.cluster.index.2 <- paste(next.level.index, 2*j, sep = "_")
      # write annotation to the next level
      iterbi.cellMeta[tmp.seuset$cellName, next.level.index] <- as.integer(tmp.seuset@active.ident) + 2*(j-1)
      # get markers
      binary_markers <- identify_binary_markers(tmp.seuset)
      b1.markers <- binary_markers[["b1"]]
      b2.markers <- binary_markers[["b2"]]
      # set end point
      if (nrow(b1.markers) >= min.marker.num) {
        b1.markers$split <- tmp.cluster.index
        b1.markers$cluster <- next.cluster.index.1
        iterbi.marker.chain <- rbind(iterbi.marker.chain, b1.markers)
      }
      if (nrow(b2.markers) >= min.marker.num) {
        b2.markers$split <- tmp.cluster.index
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



