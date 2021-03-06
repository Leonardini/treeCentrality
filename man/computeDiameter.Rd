% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Centralities.R
\name{computeDiameter}
\alias{computeDiameter}
\title{Diameter of a rooted binary phylo tree.}
\usage{
computeDiameter(tree, weight = FALSE)
}
\arguments{
\item{tree}{A phylo tree; needs to be binary and rooted, but the root can be arbitrary.}

\item{weight}{A logical scalar; if TRUE, branch lengths are used; if a root edge is present, it is ignored!}
}
\description{
\code{computeDiameter} computes the diameter of a rooted binary phylo tree.
}
\seealso{
Other network science-based tree shape statistics: \code{\link{computeBetweenness}},
  \code{\link{computeCloseness}},
  \code{\link{computeEigenvector}},
  \code{\link{computeFarness}},
  \code{\link{computeWienerIndex}}
}
\concept{network science-based tree shape statistics}
