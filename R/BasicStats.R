#' Compute the number of cherries
#'
#' \code{computeCherries} computes the number of cherries of rooted binary phylo tree
#' @param tree A phylo tree; needs to be binary and rooted.
#' @param DOUBLE A logical scalar; if TRUE, computes the number of double cherries instead.
#' @family tree shape statistics based on the number of small substructures
#' @export
computeCherries = function(tree, DOUBLE = FALSE) {
  if (is.null(tree$subtreeSizesTips)) {
    tree = addSubtreeSizes(tree)
  }
  subtreeSizes = tree$subtreeSizesTips
  goodRows = (subtreeSizes == matrix(1 + DOUBLE, nrow(subtreeSizes), ncol(subtreeSizes)))
  output = sum(rowSums(goodRows) == 2)
  output
}

#' Compute the number of pitchforks
#'
#' \code{computePitchforks} computes the number of pitchforks of rooted binary phylo tree
#' @param tree A phylo tree; needs to be binary and rooted.
#' @param FOURPRONG A logical scalar; if TRUE, computes the number of fourprongs instead.
#' @family tree shape statistics based on the number of small substructures
#' @export
computePitchforks = function(tree, FOURPRONG = FALSE) {
  if (is.null(tree$subtreeSizesTips)) {
    tree = addSubtreeSizes(tree)
  }
  subtreeS = apply(tree$subtreeSizesTips, 1, sort)
  goodCols = (subtreeS == matrix(c(1, 2 + FOURPRONG), nrow(subtreeS), ncol(subtreeS)))
  output = sum(colSums(goodCols) == 2)
  output
}

#' Compute the number of clades of sizes 4 to 8
#'
#' \code{computeNum4to8} computes the number of clades of size 4 to 8 of rooted binary phylo tree
#' @param tree A phylo tree; needs to be binary and rooted.
#' @return A vector of length 5, containing the number of clades of size 4, 5, 6, 7, and 8, in order.
#' @family tree shape statistics based on the number of small substructures
#' @export
computeNum4to8 = function(tree) {
  if (is.null(tree$subtreeSizesTips)) {
    tree = addSubtreeSizes(tree)
  }
  tab = table(rowSums(tree$subtreeSizesTips))
  output = rep(0, 5)
  allSizes = as.numeric(names(tab))
  inds = 4:8
  goodInds = which(inds %in% allSizes)
  output[inds[goodInds] - 3] =  tab[as.character(inds[goodInds])]
  output
}

#' Compute the Sackin index
#'
#' \code{computeSackin} computes the Sackin index of a rooted binary phylo tree
#' @param tree A phylo tree; needs to be binary and rooted.
#' @family tree shape statistics based on imbalance
#' @export
computeSackin = function(tree) {
  if (is.null(tree$heights)) {
    tree = addHeights(tree, weight = FALSE)
  }
  output = sum(tree$heights[1:((length(tree$heights) + 1)/2)])
  output
}

#' Compute the Colless index
#'
#' \code{computeColless} computes the Colless index of a rooted binary phylo tree
#' @param tree A phylo tree; needs to be binary and rooted.
#' @param norm A logical scalar; if TRUE, the value is normalized to lie between 0 and 1.
#' @family tree shape statistics based on imbalance
#' @export
computeColless = function(tree, norm = TRUE) {
  if (is.null(tree$subtreeSizesTips)) {
    tree = addSubtreeSizes(tree)
  }
  output = sum(abs(tree$subtreeSizesTips[,1] - tree$subtreeSizesTips[,2]))
  if (norm) {
    maxColless = choose(nrow(tree$subtreeSizesTips), 2) # attained by the caterpillar
    output = output/maxColless
  }
  output
}

#' Compute the first staircase-ness statistic
#'
#' \code{computeStairs1} computes the first staircase-ness statistic of a rooted binary phylo tree
#' @param tree A phylo tree; needs to be binary and rooted.
#' @family tree shape statistics proposed by Norstrom et al.
#' @export
computeStairs1 = function(tree) {
  if (is.null(tree$subtreeSizes)) {
    tree = addSubtreeSizes(tree)
  }
  output =  sum(tree$subtreeSizes[,1] != tree$subtreeSizes[,2])/nrow(tree$subtreeSizes)
  output
}

#' Compute the second staircase-ness statistic
#'
#' \code{computeStairs1} computes the second staircase-ness statistic of a rooted binary phylo tree
#' @param tree A phylo tree; needs to be binary and rooted.
#' @family tree shape statistics proposed by Norstrom et al.
#' @export
computeStairs2 = function(tree) {
  if (is.null(tree$subtreeSizesTips)) {
    tree = addSubtreeSizes(tree)
  }
  output = mean(apply(tree$subtreeSizesTips, 1, min)/apply(tree$subtreeSizesTips, 1, max))
  output
}

#' Compute the maximum width statistic
#'
#' \code{computeMaxWidth} computes the maximum number of nodes at a height in a rooted binary phylo tree
#' @param tree A phylo tree; needs to be binary and rooted.
#' @family tree shape statistics based on height and width
#' @export
computeMaxWidth = function(tree) {
  if (is.null(tree$heights)) {
    tree = addHeights(tree, weight = FALSE)
  }
  output = max(table(tree$heights))
  output
}

#' Compute the maximum height statistic
#'
#' \code{computeMaxHeight} computes the maximum root-to-tip distance in a rooted binary phylo tree
#' @param tree A phylo tree; needs to be binary and rooted.
#' @family tree shape statistics based on height and width
#' @export
computeMaxHeight = function(tree) {
  if (is.null(tree$heights)) {
    tree = addHeights(tree, weight = FALSE)
  }
  output = max(tree$heights)
  output
}

#' Compute the maximum delta-width statistic
#'
#' \code{computeMaxDelW} computes the maximum difference of consecutive widths in a rooted binary phylo tree
#' @param tree A phylo tree; needs to be binary and rooted.
#' @family tree shape statistics based on height and width
#' @export
computeDelW = function(tree) {
  if (is.null(tree$heights)) {
    tree = addHeights(tree, weight = FALSE)
  }
  Tab = table(tree$heights)
  output = max(diff(Tab))
  output
}
