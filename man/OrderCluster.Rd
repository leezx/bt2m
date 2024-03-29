% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.R
\name{OrderCluster}
\alias{OrderCluster}
\title{order the clusters according to similarity inside the bt2m result}
\usage{
OrderCluster(seuratObj, bt2m.result, verbose = F)
}
\arguments{
\item{seuratObj}{A Seurat object}

\item{bt2m.result}{A result file from RunBt2m() function}
}
\value{
A re-ordered list. cellMeta contains the final bifurcation for each level
marker_chain contains all the significant markers for each cluster
bifucation contains the bifurcation details (parent, child1, child2)
}
\description{
order the clusters according to similarity inside the bt2m result
}
