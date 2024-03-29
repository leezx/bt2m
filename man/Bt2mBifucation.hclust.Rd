% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bifurcated_clustering.R
\name{Bt2mBifucation.hclust}
\alias{Bt2mBifucation.hclust}
\title{Bifurcation based on hierarchical clustering}
\usage{
Bt2mBifucation.hclust(
  seuratObj,
  method = "euclidean",
  slot = "data",
  assay = "RNA"
)
}
\arguments{
\item{seuratObj}{A Seurat object}

\item{method}{The distance measurement method "euclidean or correlation"}
}
\value{
A bifurcated Seurat object (see active.ident).
}
\description{
Bifurcation based on hierarchical clustering
}
