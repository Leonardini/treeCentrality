### This function computes the height (distance to the root) of each node of a given tree
### If weight = TRUE, branch lengths are used; if a root edge is present, it is ignored!
computeHeights = function(tree, weight = FALSE) {
  stopifnot(checkPhylogeneticTree(tree))
  n = ape::Ntip(tree)
  N = 2 * n - 1
  Tab = rep(0, N)
  edges = tree$edge
  for (ind in 1:(N - 1)) {
    curRow = edges[ind,]
    Tab[curRow[2]] = ifelse(weight, tree$edge.length[ind], 1) + Tab[curRow[1]]
  }
  Tab
}

### This function computes the sizes of the left and right subtrees rooted at internal nodes
computeLRSizes = function(tree) {
  return(computeLRValues(tree, sum))
}

### This function computes the depths of the left and right subtrees rooted at internal nodes
### If weight = TRUE, the depth is the total weight of the path from the root to the node!
computeLRDepths = function(tree, weight = FALSE) {
  return(computeLRValues(tree, max, weight = weight))
}

### This function factory recursively computes values for subtrees rooted at internal nodes
computeLRValues = function(tree, FUN, weight = FALSE) {
  stopifnot(checkPhylogeneticTree(tree))
  n = ape::Ntip(tree)
  N = 2 * n - 1
  Tab = matrix(NA, n - 1, 2)
  edges = tree$edge
  for (ind in (N - 1):1) {
    curRow = edges[ind,] - n
    pos = Tab[curRow[1], 1]
    W = ifelse(weight, tree$edge.length[ind], 1)
    Tab[curRow[1], 2 - is.na(pos)] = W + ifelse(curRow[2] <= 0, 0, FUN(Tab[curRow[2],]))
  }
  Tab
}

### This function augments the tree with any subset of sizes and depths for each internal node
### and heights for each node; weight determines if the depths and heights are weighted.
augmentTree = function(tree, sizes, depths, heights, weight = FALSE) {
  if (sizes && is.null(tree$subtreeSizes)) {
    tree = addSubtreeSizes(tree)
  }
  if (depths && is.null(tree[[paste0("depths", rep("Weighted", weight))]])) {
    tree = addDepths(tree, weight)
  }
  if (heights && is.null(tree[[paste0("heights", rep("Weighted", weight))]])) {
    tree = addHeights(tree, weight)
  }
  tree
}

### This function computes the subtree sizes and adds them to a given tree
addSubtreeSizes = function(tree) {
  treeSizes = computeLRSizes(tree)
  tree$subtreeSizes = treeSizes
  tree$subtreeSizesTips  = (treeSizes + 1)/2
  tree
}

### This function computes the weighted or unweighted depths and adds them to a given tree
addDepths = function(tree, weight) {
  if (weight) {
    tree$depthsWeighted = computeLRDepths(tree, TRUE)
  } else {
    tree$depths = computeLRDepths(tree, FALSE)
  }
  tree
}

### This function computes the weighted or unweighted heights and adds them to a given tree
addHeights = function(tree, weight) {
  if (weight) {
    tree$heightsWeighted = computeHeights(tree, TRUE)
  } else {
    tree$heights = computeHeights(tree, FALSE)
  }
  tree
}

### This function checks that the input tree is phylogenetic (i.e. rooted and binary)
checkPhylogeneticTree = function(tree) {
  if (!(ape::is.binary.tree(tree)) || !(ape::is.rooted(tree))) {
    return(FALSE)
  }
  return(TRUE)
}

### This function is used to create a verbal description of a matrix based on parameters.
createName = function(Weight, Dist, Full, Lap, Norm) {
  Name = ifelse(Lap, "Laplacian", "")
  if (Dist) {
    Name = paste0("Distance", Name)
    if (Full) {
      Name = paste0(Name, ", Full")
    }
  }
  if (Lap && Norm) {
    Name = paste0(Name, ", Normalized")
  }
  if (Weight) {
    Name = paste0(Name, ", Weighted")
  }
  Name = paste(Name, "spectrum")
  Name
}
