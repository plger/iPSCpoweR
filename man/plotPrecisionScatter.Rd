% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DEA.permutations.R
\name{plotPrecisionScatter}
\alias{plotPrecisionScatter}
\title{plotPrecisionScatter}
\usage{
plotPrecisionScatter(results, plotSD = T, x = "fdr", labels = NULL,
  acol = "grey", alty = 1, alwd = 1, threshold = 0.05, ...)
}
\arguments{
\item{results}{A vector of permutation results files or the output of `readPermResults`.}

\item{plotSD}{Logical; whether to plot the standard deviation for each point (default TRUE).}

\item{x}{X axis value, either "fp" (false positives) or "fdr" (false discovery rate, default).}

\item{labels}{An optional character vector of labels to print.}

\item{acol}{color of the SD bars (default gray)}

\item{alty}{line type of the SD bars (default solid)}

\item{alwd}{line width of the SD bars (default 1)}

\item{...}{arguments passed to the plot function.}
}
\value{
Nothing, but generates a plot.
}
\description{
Plots the sensitivity and false discovery rate of a list of results of permutation analyses (i.e. permres files or output of readPermResults).
}

