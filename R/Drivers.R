basicStatNames    = c("cherries", "pitchforks", "doubcherries", "fourprong", "num4",
                      "num5", "num6", "num7", "num8", "colless", "sackin", "maxwidth",
                      "maxheight", "delW", "stairs1", "stairs2")
networkStatNames  = c("diameter", "WienerIndex", "betweenness", "closeness", "eigenvector")
spectralLMNames = c("asymmetry", "kurtosis", "densityMax", "lambdaMax")

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
#' @param unitMean A logical scalar; if TRUE, the tree's branch lengths are scaled to unit mean beforehand.
#' @return A named vector containing the 5 network science-based statistics.
#' @family drivers for computing summary statistics
#' @export
computeNetworkStats = function(tree, weight = FALSE, meanpath = FALSE, maxOnly = TRUE, unitMean = FALSE) {
  if (unitMean) {
    tree = scaleTreeBranches(tree)
  }
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
  norm = FALSE, dist = FALSE, full = FALSE, maxOnly = TRUE, unitMean = FALSE) {
  if (unitMean) {
    tree = scaleTreeBranches(tree)
  }
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

#' Summaries of the distance Lapalcian spectrum (normalized or unnormalized) of a rooted phylo tree.
#'
#' \code{computeLMStats} computes various summaries of the distance Laplacian specturm of a rooted binary phylo tree.
#' @inheritParams computeBasicStats
#' @inheritParams computeNetworkStats
#' @inheritParams computeSpectralStats
#' @return A named list containing spectral summary statistics: asymmetry, kurtosis, densityMax and lambdaMax.
#' @family drivers for computing summary statistics
#' @export
computeLMStats = function(tree, norm = FALSE, unitMean = FALSE) {
  if (unitMean) {
    tree = scaleTreeBranches(tree)
  }
  stats = spectR(tree, method = ifelse(norm, "normal", "standard"))
  output = c(stats$asymmetry, stats$peakedness1, stats$peakedness2, stats$principal_eigenvalue)
  names(output) = spectralLMNames
  output
}

#' Summaries of the distance Lapalcian spectrum (normalized or unnormalized) of a rooted phylo tree.
#'
#' \code{rankDiscriminatoryStats} ranks different summary stats of rooted binary phylo trees by discriminatory power.
#' @inheritParams computeNetworksStats
#' @param tList1 A list of trees of class \code{phylo}. The trees should be binary and rooted; else they are coerced
#' @param tList2 A list of trees of class \code{phylo}. The trees should be binary and rooted; else they are coerced
#' @param basic A logical scalar; if TRUE, computes the basic tree statistics for the lists (lengths must be equal).
#' @param basic A logical scalar; if TRUE, computes the network statistics for the lists (lengths must be equal).
#' @param spectral A logical scalar; if TRUE, computes the spectral statistics for the lists (lengths must be equal).
#' @return A named list containing the statistics for the first set of trees, the second set of trees, and p-values.
#' @family drivers for computing summary statistics
#' @export
rankDiscriminatoryStats = function(tList1, tList2, basic = TRUE, network = TRUE, spectral = TRUE, unitMean = TRUE) {
  Lens = c(length(tList1), length(tList2))
  if (Lens[1] != Lens[2]) {
    warning("Different numbers of trees are being compared")
  }
  L1 = sapply(tList1, Ntip)
  L2 = sapply(tList2, Ntip)
  if (mean(L1, na.rm = TRUE) != mean(L2, na.rm = TRUE)) {
    warning("The average numbers of tips are not equal")
  }
  nStats = length(basicStatNames) * basic + length(networkStatNames) * network + length(spectralLMNames) * spectral
  Stats = matrix(NA, Lens[1] + Lens[2], nStats)
  treeLists = c(tList1, tList2)
  for (index in 1:2) {
    startRow = ifelse(index == 1, 0, Lens[1])
    for (ind in 1:Lens[index]) {
      curCol = 0
      curTree = treeLists[[index]][[ind]]
      if (basic) {
        curRange = curCol + (1:length(basicStatNames))
        Stats[startRow + ind, curRange] = computeBasicStats(curTree, unitMean = unitMean)
        curCol = curCol + length(basicStatNames)
      }
      if (network) {
        curRange = curCol + (1:length(networkStatNames))
        Stats[startRow + ind, curRange] = computeNetworkStats(curTree, unitMean = unitMean)
        curCol = curCol + length(networkStatNames)
      }
      if (spectral) {
        curRange = curCol + (1:length(spectralLMNames))
        Stats[startRow + ind, curRange] = computeLMStats(curTree, unitMean = unitMean)
      }
    }
  }
  Stats1 = Stats[1:Lens[1], , drop = FALSE]
  Stats2 = Stats[-(1:Lens[1]), , drop = FALSE]
  comparison = compareStats(Stats1, Stats2)
  output = list(Stats1 = Stats1, Stats2 = Stats2, comparison = comparison)
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
