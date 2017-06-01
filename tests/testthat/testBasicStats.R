library(treeCentrality)
library(ape)
library(igraph)
library(phyloTop)
library(apTreeshape)
context("Basic statistics on a tree")

myTree = createTestTree()
myBasicStats = computeBasicStats(myTree)
myApTree = apTreeshape::as.treeshape(myTree)
numConfigs = phyloTop::nConfig(myTree)$numClades[as.character(4:8)]
numConfigs[is.na(numConfigs)] = 0
### Removing vector names to avoid uninformative failures
names(myBasicStats) = NULL
names(numConfigs) = NULL
uTree = makeUnweighted(myTree)
uBasicStats = computeBasicStats(uTree)
names(uBasicStats) = NULL

test_that("Basic statistics match perfectly", { # 11 tests
  expect_identical(myBasicStats, uBasicStats)
  expect_equal(myBasicStats[1], phyloTop::cherries(myTree))
  expect_equal(myBasicStats[2], phyloTop::pitchforks(myTree))
  expect_equal(myBasicStats[3] + myBasicStats[4], myBasicStats[5])
  expect_equal(myBasicStats[5:9], numConfigs)
  expect_equal(myBasicStats[10], phyloTop::colless.phylo(myTree))
  expect_equal(myBasicStats[11], apTreeshape::sackin(myApTree))
  expect_equal(myBasicStats[12], max(phyloTop::widths(myTree)))
  expect_equal(myBasicStats[13], phyloTop::maxHeight(myTree))
  expect_equal(myBasicStats[14], max(diff(phyloTop::widths(myTree))))
  expect_equal(myBasicStats[15], phyloTop::nodeImbFrac(myTree, 1))
  # expect_equal(myBasicStats[16], phyloTop::stairs(myTree)[2]) - fails; error in stairs
})
