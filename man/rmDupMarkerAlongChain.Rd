% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.R
\name{rmDupMarkerAlongChain}
\alias{rmDupMarkerAlongChain}
\title{Remove duplicated markers along a cluster chain}
\usage{
rmDupMarkerAlongChain(
  bt2m.cellMeta,
  bt2m.markerChain,
  target.cluster,
  score = "combine",
  verbose = F
)
}
\arguments{
\item{bt2m.cellMeta}{bt2m.cellMeta from bt2m}

\item{target.cluster}{target cluster}

\item{bt2m.bifucation}{bt2m.bifucation from bt2m}
}
\value{
A vector, names are cluster name, element shows the states of clusters
}
\description{
Remove duplicated markers along a cluster chain
}
