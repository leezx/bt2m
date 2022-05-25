
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

