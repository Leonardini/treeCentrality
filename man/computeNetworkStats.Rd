% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Drivers.R
\name{computeNetworkStats}
\alias{computeNetworkStats}
\title{Network science-based statistics of a rooted phylo tree.}
\usage{
computeNetworkStats(tree, weight = FALSE, meanpath = FALSE,
  maxOnly = TRUE, unitMean = FALSE)
}
\arguments{
\item{tree}{A tree of class \code{phylo}. The tree should be binary and rooted; if it is not, it will be coerced into a rooted tree using \code{root} and/or binarized using \code{multi2di}.}

\item{weight}{A logical scalar; if TRUE, the branch lengths are taken into account.}

\item{meanpath}{A logical scalar; if TRUE, the Wiener index is normalized to the mean path.}

\item{maxOnly}{A logical scalar; if TRUE, only the maximum value of each vector is returned.}

\item{unitMean}{A logical scalar; if TRUE, the tree's branch lengths are scaled to unit mean beforehand.}
}
\value{
A named vector containing the 5 network science-based statistics.
}
\description{
\code{computeNetworksStats} computes network science-based statistics of a rooted binary phylo tree.
}
\seealso{
Other drivers for computing summary statistics: \code{\link{computeBasicStats}},
  \code{\link{computeLMStats}},
  \code{\link{computeSpectralStats}},
  \code{\link{rankDiscriminatoryStats}}
}
\concept{drivers for computing summary statistics}
