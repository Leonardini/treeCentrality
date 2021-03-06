% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BasicStats.R
\name{computePitchforks}
\alias{computePitchforks}
\title{Compute the number of pitchforks}
\usage{
computePitchforks(tree, FOURPRONG = FALSE)
}
\arguments{
\item{tree}{A phylo tree; needs to be binary and rooted.}

\item{FOURPRONG}{A logical scalar; if TRUE, computes the number of fourprongs instead.}
}
\description{
\code{computePitchforks} computes the number of pitchforks of rooted binary phylo tree
}
\seealso{
Other tree shape statistics based on the number of small substructures: \code{\link{computeCherries}},
  \code{\link{computeNum4to8}}
}
\concept{tree shape statistics based on the number of small substructures}
