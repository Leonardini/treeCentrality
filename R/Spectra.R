#' Compute the matrix representation of a phylogenetic tree
#'
#' \code{computeMatrix} computes a specified matrix representation of a phylogenetic tree
#' @param tree A phylo tree; needs to be binary and rooted.
#' @param weight A logical scalar; if TRUE, branch lengths are used; if a root edge is present, it is ignored!
#' @param dist A logical scalar; if TRUE, the distance matrix is computed.
#' @param full A logical scalar; if TRUE, the node-node distances are used, otherwise, the leaf-leaf ones; ignored if dist = FALSE.
#' @param lap A logical scalar; if TRUE, the Laplacian matrix is computed (defined as diag(rowSums(M))-M for a matrix M).
#' @param norm A logical scalar; if TRUE, the Laplacian is a normalized one (i.e. the diagonal is 1); ignored if lap = FALSE.
#' @family functions used to compute tree spectra
#' @export
computeMatrix = function(tree, weight, dist, full, lap, norm) {
  stopifnot(checkPhylogeneticTree(tree))
  if (dist) {
    if (full) {
      curMatrix = ape::dist.nodes(tree)
    }
    else {
      curMatrix = ape::cophenetic.phylo(tree)
    }
  }
  else {
    Dim = 2 * tree$Nnode + 1
    curMatrix = matrix(0, Dim, Dim)
    curMatrix[tree$edge] = ifelse(rep(weight, Dim - 1), tree$edge.length, 1)
    curMatrix = curMatrix + t(curMatrix)
  }
  if (lap) {
    degrees = rowSums(curMatrix)
    curMatrix = diag(degrees) - curMatrix
    if (norm) {
      curMatrix = curMatrix/sqrt(outer(degrees, degrees))
    }
  }
  curMatrix
}

#' Compute the spectrum (eigenvalues) of a matrix representing a phylogenetic tree
#'
#' \code{computeSpectrum} computes the spectrum of a specified matrix of a phylogenetic tree
#' @param tree A phylo tree; needs to be binary and rooted.
#' @param weight A logical scalar; if TRUE, branch lengths are used; if a root edge is present, it is ignored!
#' @param dist A logical scalar; if TRUE, the distance matrix is computed.
#' @param full A logical scalar; if TRUE, the node-node distances are used, otherwise, the leaf-leaf ones; ignored if dist = FALSE.
#' @param lap A logical scalar; if TRUE, the Laplacian matrix is computed (defined as diag(rowSums(M))-M for a matrix M).
#' @param norm A logical scalar; if TRUE, the Laplacian is a normalized one (i.e. the diagonal is 1); ignored if lap = FALSE.
#' @family functions used to compute tree spectra
#' @export
computeSpectrum  = function(tree, weight, dist, full, lap, norm) {
  curMatrix = computeMatrix(tree, weight, dist, full, lap, norm)
  curSpectrum = eigen(curMatrix, symmetric = TRUE, only.values = TRUE)$values
  curSpectrum
}
