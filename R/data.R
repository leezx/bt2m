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