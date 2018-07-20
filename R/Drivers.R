basicStatNames    = c("cherries", "pitchforks", "doubcherries", "fourprong", "num4",
                      "num5", "num6", "num7", "num8", "colless", "sackin", "maxwidth",
                      "maxheight", "delW", "stairs1", "stairs2")
networkStatNames  = c("diameter", "WienerIndex", "betweenness", "closeness", "eigenvector")

#' Basic summary statistics of a rooted phylo tree.
#'
#' \code{computeBasicStats} computes basic summary statistics of a rooted binary phylo tree.
#' @param tree A tree of class \code{phylo}. The tree should be binary and rooted; if it is not, it will be coerced into a rooted tree using \code{root} and/or binarized using \code{multi2di}.
#' @return A named vector containing the 16 basic summary statistics.
#' @family drivers for computing summary statistics
#' @export
computeBasicStats = function(tree) {
  tree = augmentTree(tree, sizes = TRUE, depths = FALSE, heights = TRUE, weight = FALSE)
  cherries      = computeCherries(tree, FALSE)
  doubcherries  = computeCherries(tree, TRUE)
  pitchforks    = computePitchforks(tree, FALSE)
  fourprong     = computePitchforks(tree, TRUE)
  num4to8       = computeNum4to8(tree)
  colless       = computeColless(tree)
  sackin        = computeSackin(tree)
  maxWidth      = computeMaxWidth(tree)
  maxHeight     = computeMaxHeight(tree)
  delW          = computeDelW(tree)
  stairs1       = computeStairs1(tree)
  stairs2       = computeStairs2(tree)
  output = c(cherries, pitchforks, doubcherries, fourprong, num4to8,
             colless, sackin, maxWidth, maxHeight, delW, stairs1, stairs2)
  names(output) = basicStatNames
  output
}

#' Network science-based statistics of a rooted phylo tree.
#'
#' \code{computeNetworksStats} computes network science-based statistics of a rooted binary phylo tree.
#' @inheritParams computeBasicStats
#' @param weight A logical scalar; if TRUE, the branch lengths are taken into account.
#' @param meanpath A logical scalar; if TRUE, the Wiener index is normalized to the mean path.
#' @param maxOnly A logical scalar; if TRUE, only the maximum value of each vector is returned.
#' @return A named vector containing the 5 network science-based statistics.
#' @family drivers for computing summary statistics
#' @export
computeNetworkStats = function(tree, weight = FALSE, meanpath = FALSE, maxOnly = TRUE) {
  tree = augmentTree(tree, sizes = TRUE, depths = TRUE, heights = TRUE, weight = weight)
  diameter      = computeDiameter(tree, weight = weight)
  WienerIndex   = computeWienerIndex(tree, norm = meanpath, weight = weight)
  betweenness   = computeBetweenness(tree, weight = weight)
  closeness     = computeCloseness(tree, weight = weight)
  eigenvector   = computeEigenvector(tree, weight = weight, scale = FALSE)[[1]]
  if (maxOnly) {
    betweenness = max(betweenness)
    closeness   = max(closeness)
    eigenvector = max(eigenvector)
  }
  output = c(diameter, WienerIndex, betweenness, closeness, eigenvector)
  L = length(eigenvector)
  names(output) = c(networkStatNames[1:2], rep(networkStatNames[-(1:2)], each = L))
  output
}

#' Spectra of a rooted phylo tree.
#'
#' \code{computeSpectralStats} computes various spectra of a rooted binary phylo tree.
#' @inheritParams computeBasicStats
#' @inheritParams computeNetworkStats
#' @param weight A logical vector; a TRUE (FALSE) entry means the branch lengths are (not) taken into account.
#' @param adj A logical vector; a TRUE (FALSE) entry means the adjacency (Laplacian) spectrum is computed.
#' @param norm A logical vector; a TRUE (FALSE) entry means the Laplacian is (not) normalized.
#' @param dist A logical vector; a TRUE (FALSE) entry means the distance (regular) matrix is used.
#' @param full A logical vector;  a TRUE (FALSE) entry means the full (leaf-leaf) distance matrix is used.
#' @return A named list containing the specified spectra (or only their maximum values).
#' @family drivers for computing summary statistics
#' @export
computeSpectralStats = function(tree, weight = c(FALSE, TRUE), adj = c(FALSE, TRUE),
  norm = FALSE, dist = FALSE, full = FALSE, maxOnly = TRUE) {
  tree = checkPhylogeneticTree(tree)
  normFactor = sum(ifelse(adj, 1, length(norm))) # norm only applies to adj = FALSE (lap)
  distFactor = sum(ifelse(dist, length(full), 1)) # full only applied to dist = TRUE
  numSpectra = length(weight) * length(adj) * normFactor * distFactor
  allSpectra = vector("list", numSpectra)
  ind = 0
  for (Weight in weight) {
    for (Adj in adj) {
      for (Norm in norm) {
        if (Adj && Norm) { # validity check: Norm only makes sense if Adj = FALSE
          next
        }
        for (Dist in dist) {
          for (Full in full) {
            if (!Dist && Full) {  # validity check: Full only makes sense if Dist = TRUE
              next
            }
            curName = createName(Weight, Dist, Full, !Adj, Norm)
            allSpectra[[curName]] = computeSpectrum(tree, Weight, Dist, Full, !Adj, Norm)
            ind = ind + 1
          }
        }
      }
    }
  }
  allSpectra = allSpectra[!sapply(allSpectra,is.null)]
  if (maxOnly) {
    allSpectra = sapply(allSpectra, max)
  }
  allSpectra
}

### TODO: ADD DESCRIPTION!
rankDiscriminatoryStats = function(treeList1, treeList2, basic = TRUE, network = TRUE, spectral = TRUE) {
  Lens = c(length(treeList1), length(treeList2))
  if (Lens[1] != Lens[2]) {
    warning("Different numbers of trees are being compared")
  }
  L1 = sapply(treeList1, Ntip)
  L2 = sapply(treeList2, Ntip)
  if (mean(L1, na.rm = TRUE) != mean(L2, na.rm = TRUE)) {
    warning("The average numbers of tips are not equal")
  }
  numStats = length(basicStatNames) * basic + length(networkStatNames) * network + 4 * spectral
  Stats = matrix(NA, Lens[1] + Lens[2], numStats)
  for (index in 1:2) {
    startRow = ifelse(index == 1, 0, Lens[1])
    for (ind in 1:Lens[index]) {
      curCol = 0
      curTree = treeLists[[index]][[ind]]
      if (basic) {
        Stats[startRow + ind, curCol + (1:length(basicStatNames))] = computeBasicStats(curTree)
        curCol = curCol + length(basicStatNames)
      }
      if (network) {
        Stats[startRow + ind, curCol + (1:length(networkStatNames))] = computeNetworkStats(curTree)
        curCol = curCol + length(networkStatNames)
      }
      if (spectral) {
        
      }
    }
  }
  Stats1 = Stats[1:Lens[1], , drop = FALSE]
  Stats2 = Stats[-(1:Lens[1]), , drop = FALSE]
  output = compareStats(Stats1, Stats2)
  output
}

### TODO: ADD DESCRIPTION!
compareStats = function(Stats1, Stats2) {
  stopifnot(ncol(Stats1) == ncol(Stats2))
  stopifnot(all(sort(colnames(Stats1)) == sort(colnames(Stats2))))
  Mat = matrix(NA, ncol(Stats1), 3, dimnames = list(colnames(Stats1), c("mean1", "mean2", "p-value")))
  for (ind in 1:ncol(Stats1)) {
    curName = colnames(Stats1)[ind]
    cur1 = Stats1[,curName]
    cur2 = Stats2[,curName]
    curTest = t.test(unlist(cur1), unlist(cur2))
    curResult = c(curTest$estimate, curTest$p.value)
    names(curResult) = NULL
    Mat[ind,] = curResult
  }
  Mat = Mat[order(Mat[, "p-value"]), ]
  Mat
}