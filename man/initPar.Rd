% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getFunctions.R
\name{initPar}
\alias{initPar}
\title{initPar}
\usage{
initPar(ncores = NULL)
}
\arguments{
\item{ncores}{The number of cores to use, default detectedCores()-1.}
}
\value{
A makeCluster object if parallel is supported and more than one core is used; otherwise FALSE.
}
\description{
Initializes multithreading. This does not need to be called directly; the package's functions will call this.
}
