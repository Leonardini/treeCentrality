#' Diameter of a rooted binary phylo tree.
#'
#' \code{computeDiameter} computes the diameter of a rooted binary phylo tree.
#'
#' @param tree A phylo tree; needs to be binary and rooted, but the root can be arbitrary.
#' @param weight A logical scalar; if TRUE, branch lengths are used; if a root edge is present, it is ignored!
#' @family network science-based tree shape statistics
#' @export
computeDiameter = function(tree, weight = FALSE) {
  if (is.null(tree[[paste0("depths", rep("Weighted", weight))]])) {
    tree = addDepths(tree, weight = weight)
  }
  diam = max(rowSums(tree[[paste0("depths", rep("Weighted", weight))]]))
  diam
}

#' Wiener index of a rooted binary phylo tree.
#'
#' \code{computeWienerIndex} computes the Wiener index of a rooted binary phylo tree.
#'
#' @param tree A phylo tree; needs to be binary and rooted, but the root can be arbitrary.
#' @param norm A logical scalar; if TRUE, the Wiener index is normalized to give the average path.
#' @param weight A logical scalar; if TRUE, branch lengths are used; if a root edge is present, it is ignored!
#' @family network science-based tree shape statistics
#' @export
computeWienerIndex = function(tree, norm = FALSE, weight = FALSE) {
  if (is.null(tree$subtreeSizes)) {
    tree = addSubtreeSizes(tree)
  }
  q = rowSums(tree$subtreeSizes) + 1
  n = length(q)
  N = 2 * n + 1
  stopifnot(q[1] == N)
  if (weight) {
    edges = tree$edge
    lengths = tree$edge.length
    stopifnot(all(lengths >= 0))
    W = 0
    for (ind in 1:nrow(tree$edge)) {
      curEndpoint = edges[ind, 2]
      curQ = ifelse(curEndpoint > (n + 2), q[curEndpoint - (n + 1)], 1)
      W = W + curQ * (N - curQ) * lengths[ind]
    }
  }
  else {
    W = sum(q * (N - q)) + (N - 1) * (n + 1)
  }
  W = W / ifelse(norm, choose(N, 2), 1)
  W
}

#' Betweenness centrality of a rooted binary phylo tree.
#'
#' \code{computeBetweenness} computes the betweenness centrality of a rooted binary phylo tree.
#' @param tree A phylo tree; needs to be binary and rooted, but the root can be arbitrary.
#' @param weight A logical scalar; if TRUE, branch lengths are used; however, they do not change the result!
#' @family network science-based tree shape statistics
#' @export
computeBetweenness = function(tree, weight = FALSE) {
  if (is.null(tree$subtreeSizes)) {
    tree = addSubtreeSizes(tree)
  }
  Tab = tree$subtreeSizes
  n = nrow(Tab)
  rSums = rowSums(Tab)
  Centralities = c(rep(0, n + 1), Tab[,1] * Tab[,2] + rSums * (2 * n - rSums))
  Centralities
}

#' Closeness centrality of a rooted binary phylo tree.
#'
#' \code{computeCloseness} computes the closeness centrality of a rooted binary phylo tree.
#' @param tree A phylo tree; needs to be binary and rooted, but the root can be arbitrary.
#' @param weight A logical scalar; if TRUE, branch lengths are used; if a root edge is present, it is ignored!
#' @family network science-based tree shape statistics
#' @export
computeCloseness = function(tree, weight = FALSE) {
  return(1/computeFarness(tree, weight = weight))
}

#' Farness of a rooted binary phylo tree.
#'
#' \code{computeFarness} computes the faress of each node in a rooted binary phylo tree.
#' @param tree A phylo tree; needs to be binary and rooted, but the root can be arbitrary.
#' @param weight A logical scalar; if TRUE, branch lengths are used; if a root edge is present, it is ignored!
#' @family network science-based tree shape statistics
#' @export
computeFarness = function(tree, weight = FALSE) {
  if (is.null(tree$subtreeSizes)) {
    tree = addSubtreeSizes(tree)
  }
  sizes = rowSums(tree$subtreeSizes)
  n = ape::Ntip(tree)
  N = 2 * n - 1
  Farness = rep(NA, N)
  if (weight) {
    if (is.null(tree$heightsWeighted)) {
      tree = addHeights(tree, weight = TRUE)
    }
    Farness[n + 1] = sum(tree$heightsWeighted) ### bug fix: originally had only sum[1:n]
  } else {
    Farness[n + 1] = sum(sizes)
  }
  edges = tree$edge
  for (ind in 1:nrow(edges)) {
    curRow = edges[ind,]
    kid = curRow[2]
    subSize = 1 + ifelse(kid <= n, 0, sizes[kid - n])
    W = ifelse(weight, tree$edge.length[ind], 1)
    Farness[kid] = Farness[curRow[1]] + (N - 2 * subSize) * W
  }
  Farness
}

#' Eigenvector centrality of a rooted binary phylo tree.
#'
#' \code{computeEigenvector} computes the eigenvector centrality of a rooted binary phylo tree.
#' @param tree A phylo tree; needs to be binary and rooted, but the root can be arbitrary.
#' @param weight A logical scalar; if TRUE, branch lengths are used; if a root edge is present, it is ignored!
#' @param scale A logical scalar; if TRUE, the resulting vector is scaled to have a maximum value of 1.
#' If no scaling is used (the default), the eigenvector has unit length in the Euclidean norm.
#' @return A list containing the dominant eigenvector and its corresponding eigenvalue.
#' @family network science-based tree shape statistics
#' @export
computeEigenvector = function(tree, weight = FALSE, scale = FALSE) {
  stopifnot(checkPhylogeneticTree(tree))
  graph = ape::as.igraph.phylo(tree, directed = FALSE)
  igraph::E(graph)$weight = ifelse(rep(weight, nrow(tree$edge)), tree$edge.length, 1)
  adj_matrix = igraph::get.adjacency(graph, sparse = TRUE, attr = "weight")
  EV = rARPACK::eigs_sym(adj_matrix, k = 1, which = "LM", opts = list(retvec = TRUE))
  evector = abs(EV$vectors[,1])
  if (scale) {
    evector = evalue/max(evector)
  }
  names(evector) = rownames(adj_matrix)
  evalue = abs(EV$values[1])
  orderedVector = rep(NA, length(evector))
  n = ape::Ntip(tree)
  orderedVector[1:n] = evector[tree$tip.label]
  nodeLabels = tree$node.label
  nodeNames = ifelse(rep(is.null(nodeLabels), n-1), paste0("Node", 1:(n-1)), nodeLabels)
  orderedVector[-(1:n)] = evector[nodeNames]
  output = list(eigenvector = orderedVector, eigenvalue = evalue)
  output
}
