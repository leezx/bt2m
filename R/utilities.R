
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

add.missing.DEGs <- function(epxrM, all.genes) {
  # to matrix
  #epxrM <- as.matrix(epxrM)
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
  #epxrM <- as(epxrM, "dgCMatrix")
  return(epxrM)
}

#' a subset function for seurat obj, origninal cannot work sometimes
subset_cells <- function(seuratObj, condition) {
    names(condition) <- colnames(seuratObj)
    seuratObj$select_cells <- factor(condition, levels = c(T, F))
    tmp.seuratObj <- subset(seuratObj, subset = select_cells == T)
    return(tmp.seuratObj)
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

#' A general function to identify DEGs between case and control
DEG.1by1 <- function(seuratObj, ident.1 = "Vcl cKO", ident.2 = "Control", assay = "RNA") {
    DEGs <- list()
    for (i in unique(seuratObj@active.ident)) {
        message(sprintf("identifying DEGs in %s between %s and %s...", i, ident.1, ident.2))
        tmp.seuratObj <- subset(seuratObj, subset = cluster == i)
        tmp.seuratObj@active.ident <- tmp.seuratObj$group
        tmp.DEGs <- FindMarkers(tmp.seuratObj, ident.1 = ident.1, ident.2 = ident.2, only.pos = F, 
                                min.pct = 0, min.diff.pct = "-Inf", logfc.threshold = 0, assay = assay)
        if(dim(tmp.DEGs)[1]<1) {
          message(sprintf("skipping %s, no DEGs found...", i))
          next
        }
        DEGs[[i]] <- add.missing.DEGs(tmp.DEGs, rownames(tmp.seuratObj@assays$RNA@counts))
        # fill emplty genes
        tmp.empty <- rowSums(DEGs[[i]])==0
        DEGs[[i]][tmp.empty,]$p_val <- 1
        DEGs[[i]][tmp.empty,]$p_val_adj <- 1
        # add gene column
        DEGs[[i]]$gene <- rownames(DEGs[[i]])
        # add log2FC and correlation
        cells_1 <- rownames(subset(seuratObj@meta.data, cluster==i & group==ident.2))
        cells_2 <- rownames(subset(seuratObj@meta.data, cluster==i & group==ident.1))
        log2fc <- log2(rowMeans(seuratObj@assays$RNA@data[,cells_2])+1) - log2(rowMeans(seuratObj@assays$RNA@data[,cells_1])+1)
        corM <- cor(t(as.matrix(seuratObj@assays$RNA@data[,c(cells_1,cells_2)])), c(rep(0, length.out = length(cells_1)),
                                                                       rep(1, length.out = length(cells_2))
                                                                       ))
        corM[is.na(corM)] <- 0
        DEGs[[i]]$log2FC <- log2fc[DEGs[[i]]$gene]
        DEGs[[i]]$correlation <- as.data.frame(corM)[DEGs[[i]]$gene,]
        #
    }
    return(DEGs)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions for plotting
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' origninal DotPlot in Seurat cannot order features

source("https://github.com/satijalab/seurat/raw/8da35ba8dc2d60f504f7cf276efd684dc19a418c/R/utilities.R")
DotPlot_order <- function (object, assay = NULL, features, cols = c("lightgrey", 
    "blue"), col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 6, 
    group.by = NULL, split.by = NULL, scale = TRUE, scale.by = "radius", 
    scale.min = NA, scale.max = NA) 
{
    if(is.null(assay)) assay <- DefaultAssay(object = object)
    DefaultAssay(object = object) <- assay
    scale.func <- switch(EXPR = scale.by, size = scale_size, 
        radius = scale_radius, stop("'scale.by' must be either 'size' or 'radius'"))
    data.features <- FetchData(object = object, vars = features)
    data.features$id <- if (is.null(x = group.by)) {
        Idents(object = object)
    }
    else {
        object[[group.by, drop = TRUE]]
    }
    if (!is.factor(x = data.features$id)) {
        data.features$id <- factor(x = data.features$id)
    }
    id.levels <- levels(x = data.features$id)
    data.features$id <- as.vector(x = data.features$id)
    if (!is.null(x = split.by)) {
        splits <- object[[split.by, drop = TRUE]]
        if (length(x = unique(x = splits)) > length(x = cols)) {
            stop("Not enought colors for the number of groups")
        }
        cols <- cols[1:length(x = unique(x = splits))]
        names(x = cols) <- unique(x = splits)
        data.features$id <- paste(data.features$id, splits, sep = "_")
        unique.splits <- unique(x = splits)
        id.levels <- paste0(rep(x = id.levels, each = length(x = unique.splits)), 
            "_", rep(x = unique(x = splits), times = length(x = id.levels)))
    }
    data.plot <- lapply(X = unique(x = data.features$id), FUN = function(ident) {
        data.use <- data.features[data.features$id == ident, 
            1:(ncol(x = data.features) - 1), drop = FALSE]
        avg.exp <- apply(X = data.use, MARGIN = 2, FUN = function(x) {
            return(mean(x = expm1(x = x)))
        })
        pct.exp <- apply(X = data.use, MARGIN = 2, FUN = PercentAbove, 
            threshold = 0)
        return(list(avg.exp = avg.exp, pct.exp = pct.exp))
    })
    names(x = data.plot) <- unique(x = data.features$id)
    data.plot <- lapply(X = names(x = data.plot), FUN = function(x) {
        data.use <- as.data.frame(x = data.plot[[x]])
        data.use$features.plot <- rownames(x = data.use)
        data.use$id <- x
        return(data.use)
    })
    data.plot <- do.call(what = "rbind", args = data.plot)
    if (!is.null(x = id.levels)) {
        data.plot$id <- factor(x = data.plot$id, levels = id.levels)
    }
    if (length(x = levels(x = data.plot$id)) == 1) {
        scale <- FALSE
        warning("Only one identity present, the expression values will be not scaled.")
    }
    avg.exp.scaled <- sapply(X = unique(x = data.plot$features.plot), 
        FUN = function(x) {
            data.use <- data.plot[data.plot$features.plot == 
                x, "avg.exp"]
            if (scale) {
                data.use <- scale(x = data.use)
                data.use <- MinMax(data = data.use, min = col.min, 
                  max = col.max)
            }
            else {
                data.use <- log(x = data.use)
            }
            return(data.use)
        })
    avg.exp.scaled <- as.vector(x = t(x = avg.exp.scaled))
    if (!is.null(x = split.by)) {
        avg.exp.scaled <- as.numeric(x = cut(x = avg.exp.scaled, 
            breaks = 20))
    }
    data.plot$avg.exp.scaled <- avg.exp.scaled
    data.plot$features.plot <- factor(x = data.plot$features.plot, 
        levels = rev(x = features))
    data.plot$pct.exp[data.plot$pct.exp < dot.min] <- NA
    data.plot$pct.exp <- data.plot$pct.exp * 100
    if (!is.null(x = split.by)) {
        splits.use <- vapply(X = as.character(x = data.plot$id), 
            FUN = gsub, FUN.VALUE = character(length = 1L), pattern = paste0("^((", 
                paste(sort(x = levels(x = object), decreasing = TRUE), 
                  collapse = "|"), ")_)"), replacement = "", 
            USE.NAMES = FALSE)
        data.plot$colors <- mapply(FUN = function(color, value) {
            return(colorRampPalette(colors = c("grey", color))(20)[value])
        }, color = cols[splits.use], value = avg.exp.scaled)
    }
    color.by <- ifelse(test = is.null(x = split.by), yes = "avg.exp.scaled", 
        no = "colors")
    if (!is.na(x = scale.min)) {
        data.plot[data.plot$pct.exp < scale.min, "pct.exp"] <- scale.min
    }
    if (!is.na(x = scale.max)) {
        data.plot[data.plot$pct.exp > scale.max, "pct.exp"] <- scale.max
    }
    #
    # print(head(data.plot))
    data.plot$features.plot <- factor(data.plot$features.plot, levels = features)
    #
    plot <- ggplot(data = data.plot, mapping = aes_string(x = "features.plot", 
        y = "id")) + geom_point(mapping = aes_string(size = "pct.exp", 
        color = color.by)) + scale.func(range = c(0, dot.scale), 
        limits = c(scale.min, scale.max)) + theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank()) + guides(size = guide_legend(title = "Percent Expressed")) + 
        labs(x = "Features", y = ifelse(test = is.null(x = split.by), 
            yes = "Identity", no = "Split Identity")) + theme_cowplot()
    if (!is.null(x = split.by)) {
        plot <- plot + scale_color_identity()
    }
    else if (length(x = cols) == 1) {
        plot <- plot + scale_color_distiller(palette = cols)
    }
    else {
        plot <- plot + scale_color_gradient(low = cols[1], high = cols[2])
    }
    if (is.null(x = split.by)) {
        plot <- plot + guides(color = guide_colorbar(title = "Average Expression"))
    }
    return(plot)
}


################################################################
#' a stacked version of VlnPlot in Seurat, from public website

## remove the x-axis text and tick
## plot.margin to adjust the white space between each plot.
## ... pass any arguments to VlnPlot in Seurat
modify_vlnplot <- function(obj, 
                          feature, 
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  p <- VlnPlot(obj, features = feature, pt.size = pt.size, ... )  + 
                xlab("") + ylab(feature) + ggtitle("") + 
                theme(legend.position = "none", 
                      axis.text.x = element_blank(), 
                      axis.ticks.x = element_blank(), 
                      axis.title.y = element_text(size = 20, angle = 0, face = "bold.italic", vjust = 0.5), 
                      axis.text.y = element_text(size = 10), 
                      plot.margin = plot.margin ) 
  return(p)
}

## extract the max value of the y axis
extract_max <- function(p){
  ymax <- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}


## main function
StackedVlnPlot <- function(obj, features,
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  
  plot_list <- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
  
  # Add back x-axis title to bottom plot. patchwork is going to support this?
  plot_list[[length(plot_list)]] <- plot_list[[length(plot_list)]] +
    theme(axis.text.x = element_text(angle = 75, size = 15), axis.ticks.x = element_line())
  
  # change the y-axis tick to only max value 
  ymaxs <- purrr::map_dbl(plot_list, extract_max)
  plot_list <- purrr::map2(plot_list, ymaxs, function(x,y) x + 
                            scale_y_continuous(breaks = c(y)) + 
                            expand_limits(y = y))

  p <- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  # p <- cowplot::plot_grid(plotlist = plot_list, ncol = 1, labels = features) # not aligned
  return(p)
}

################################################################


