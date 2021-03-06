% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DEA.permutations.R
\name{getSensitivityMatrices}
\alias{getSensitivityMatrices}
\title{getSensitivityMatrices}
\usage{
getSensitivityMatrices(reslist, bins = 5, unlogExpr = T,
  display_numbers = T)
}
\arguments{
\item{reslist}{A list of analysis results, such as produced by the `readPermResults()` function.}

\item{bins}{Either an integer indicating the number of bins in which to (try to) split the expression levels (default 5), or a list of ranges for binning.}

\item{unlogExpr}{Whether to "un-log" read counts for labels (default TRUE).}
}
\value{
a heatmap.
}
\description{
Plot a joined sensitivity matrix across foldchanges and expression levels
}

