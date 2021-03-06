% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DEA.permutations.R
\name{voomWrapperSumReps}
\alias{voomWrapperSumReps}
\title{voomWrapperSumReps}
\usage{
voomWrapperSumReps(e, mm = NULL, nested = NULL, DEsamples = NULL,
  paired = NULL)
}
\arguments{
\item{e}{the expression matrix}

\item{mm}{the design model (as created by the `model.matrix()' function), for paired analyses only. Default NULL.}

\item{nested}{NULL or character vector of size=ncol(e), indicating the blocking variable (default NULL).}

\item{DEsamples}{the index of the samples in which the differential expression was added. Mandatory if mm=NULL, ignored otherwise.}

\item{paired}{Ignored; for consistency with edgeRwrapper.}
}
\value{
A data.frame.
}
\description{
Performs one differential expression analysis using voom
}

