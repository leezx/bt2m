% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.R
\name{RunBT2M.DC.extend}
\alias{RunBT2M.DC.extend}
\title{The main function to perform iteratively bifurcation clustering RunBT2M.DC.extend}
\usage{
RunBT2M.DC.extend(
  seuratObj,
  bt2m.result,
  slot = "data",
  assay = "RNA",
  method = "graph",
  min.marker.num = 100,
  max.level.num = 100,
  min.cell.count = 50,
  verbose = F
)
}
\arguments{
\item{seuratObj}{A Seurat object}

\item{method}{The method to perform bifurcation clustering "graph (default), hclust or kmeans"}

\item{min.marker.num}{Minimal number of markers to confirm a bifurcation}

\item{max.level.num}{Maximum number of level for bifurcation}

\item{min.cell.count}{Minimal number of cells to perform bifurcation (must bigger than PC number: 50)}

\item{verbose}{Print detail proccessing messages}

\item{resolution.sets}{The number of resolution for searching}
}
\value{
A list. cellMeta contains the preliminary bifurcation for each level
marker_chain contains all the significant markers for each cluster
bifucation contains the bifurcation details (parent, child1, child2)
}
\description{
The main function to perform iteratively bifurcation clustering RunBT2M.DC.extend
}
