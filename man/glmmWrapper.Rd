% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wrappers.R
\name{glmmWrapper}
\alias{glmmWrapper}
\title{glmmWrapper}
\usage{
glmmWrapper(e, nested, DEsamples = NULL, mm = NULL, paired = NULL)
}
\arguments{
\item{e}{the expression matrix}

\item{nested}{character vector of size=ncol(e), indicating the blocking variable}

\item{DEsamples}{the index of the samples in which the differential expression was added. Mandatory if mm=NULL, ignored otherwise.}

\item{mm}{Ignored; for consistency with alternative functions.}

\item{paired}{Ignored; for consistency with alternative functions.}
}
\value{
A data.frame.
}
\description{
Performs one differential expression analysis using mixed models (via `MASS::glmmPQL`) and quasipoisson distribution.
}
