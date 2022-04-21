#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Download raw single-cell matrix for testing
#'
#' @return dowloaded single-cell dataset
#' @export
#'
GetTestData <- function() {
    # Seurat test data
    download.file("https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz",
                  destfile="pbmc3k_filtered_gene_bc_matrices.tar.gz")
    untar("pbmc3k_filtered_gene_bc_matrices.tar.gz")
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Data
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Human transcription factors
#'
#' A vector of genes used in marker annotation
#'
#' @format A vector
#' @concept data
#' @source \url{http://bioinfo.life.hust.edu.cn/HumanTFDB/#!/download}
#'
"human.tfs"

#' Mouse transcription factors
#'
#' A vector of genes used in marker annotation
#'
#' @format A vector
#' @concept data
#' @source \url{http://bioinfo.life.hust.edu.cn/AnimalTFDB/#!/download}
#'
"mouse.tfs"

#' Human ligand and receptor
#'
#' A list of genes used in marker annotation
#'
#' @format A list of two vectors
#' \describe{
#'   \item{human.ligand}{ligand genes in human}
#'   \item{human.receptor}{receptor genes in human}
#' }
#' @concept data
#' @source \url{http://www.cellchat.org/}
#'
"human.ligand.receptor"

#' Mouse ligand and receptor
#'
#' A list of genes used in marker annotation
#'
#' @format A list of two vectors
#' \describe{
#'   \item{mouse.ligand}{ligand genes in mouse}
#'   \item{mouse.receptor}{receptor genes in mouse}
#' }
#' @concept data
#' @source \url{http://www.cellchat.org/}
#'
"mouse.ligand.receptor"
