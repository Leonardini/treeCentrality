### This function computes the height (distance to the root) of each node of a given tree
### If weight = TRUE, branch lengths are used; if a root edge is present, it is ignored!
computeHeights = function(tree, weight = FALSE) {
  tree = checkPhylogeneticTree(tree)
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
  tree = checkPhylogeneticTree(tree)
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
### If not, arbitrarily roots it at the 1st internal node and/or binarizes it by multi2di
checkPhylogeneticTree = function(tree) {
  if (!(ape::is.rooted(tree))) {
    tree = ape::root(tree, node = ape::Ntip(tree) + 2, resolve.root = TRUE)
  }
  if (!(ape::is.binary.tree(tree))) {
    tree = ape::multi2di(tree)
  }
  return(tree)
}

### This function is used to create a verbal description of a matrix based on parameters.
createName = function(Weight, Dist, Full, Lap, Norm) {
  Name = ifelse(Lap, "Laplacian", "")
  if (Dist) {
    Name = paste0("Distance", ifelse(Name == "", "", " "), Name)
    if (Full) {
      Name = paste0(Name, ", Full")
    }
  }
  if (Lap && Norm) {
    Name = paste0(Name, ", Normalized")
  }
  if (Weight) {
    Name = paste0(Name, ifelse(Name == "", "", ", "), "Weighted")
  }
  Name = paste(Name, "spectrum")
  Name
}

### Utility function from the NeighborJoin.R file
createEmptyTree = function(numNodes) {
  Tree = rep(list(c(NA, Inf, NA, NA)), numNodes)
  names(Tree) = 1:numNodes
  Tree
}

### Utility function from the NeighborJoin.R file
createNewNode = function(Tree, children, branchLengths, position) {
  N = position
  for (ind in 1:length(children)) {
    curChild = children[ind]
    Tree[[curChild]][1] = N
    Tree[[curChild]][2] = branchLengths[ind]
  }
  Tree[[N]] = c(NA, Inf, children)
  Tree
}

### Conversion function from the NeighborJoin.R file
convertPhylo = function(Tree) {
  N = length(Tree)
  n = (N + 1)/2
  Labels = names(Tree)[1:n]
  Edges = matrix(NA, N - 1, 3)
  for (ind in 1:(N - 1)) {
    curElement = Tree[[ind]]
    Edges[ind,] = c(curElement[1], ind, curElement[2])
  }
  Edges[Edges[,1] > n, 1] = 3 * n - Edges[Edges[,1] > n, 1]
  Edges[Edges[,2] > n, 2] = 3 * n - Edges[Edges[,2] > n, 2]
  iTree = igraph::graph(t(Edges[, 1:2]), directed = TRUE)
  orderT = igraph::graph.dfs(iTree, root = 3 * n - N)$order
  Edges = Edges[order(Edges[,2]),]
  newEdges = matrix(NA, nrow(Edges), ncol(Edges))
  newEdges[order(orderT[-1]),] = Edges
  newLengths = newEdges[,3]
  newEdges = newEdges[,-3]
  randTree = ape::rtree(n, rooted = TRUE)
  randTree$edge = newEdges
  randTree$edge.length = newLengths
  randTree$tip.label = Labels
  randTree
}

### This function creates the tree used to cross-check the statistics
createTestTree = function() {
  Tr = createEmptyTree(numNodes = 13)
  Tr = createNewNode(Tr, c(4,5), c(1,1), 8)
  Tr = createNewNode(Tr, c(6,7), c(1,1), 9)
  Tr = createNewNode(Tr, c(2,3),  c(2,2), 10)
  Tr = createNewNode(Tr, c(8,9), c(1,1), 11)
  Tr = createNewNode(Tr, c(10,11), c(2,2), 12)
  Tr = createNewNode(Tr, c(1,12), c(8,4), 13)
  Tree = convertPhylo(Tr)
  Tree$tip.label = c("a","f","g","d","e","b","c")
  Tree$node.label = LETTERS[1:6]
  Tree
}

### This function creates an unweighted version of a weighted tree
makeUnweighted = function(Tree) {
  uTree = Tree
  uTree$edge.length = rep(1, nrow(Tree$edge))
  uTree
}
