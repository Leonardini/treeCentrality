library(treeCentrality)
library(ape)
library(igraph)
context("Network statistics on a tree")

myTree = createTestTree()
myNetStats = computeNetworkStats(myTree, weight = TRUE, meanpath = TRUE, maxOnly = FALSE, unitMean = FALSE)
myDF = as.data.frame(cbind(myTree$edge, weight = myTree$edge.length))
myG = igraph::graph_from_data_frame(myDF, directed = FALSE)
N = igraph::gorder(myG)
myDistances = distances(myG)
diag(myDistances) = NA
myBetweenness = igraph::betweenness(myG, directed = FALSE)
order = order(as.numeric(names(myBetweenness)))
myBetweenness = myBetweenness[order]
myCloseness = igraph::closeness(myG)
myCloseness = myCloseness[order]
myEigenvector = igraph::eigen_centrality(myG, scale = FALSE)$vector
myEigenvector = myEigenvector[order]
names(myNetStats) = NULL
names(myBetweenness) = NULL
names(myCloseness) = NULL
names(myEigenvector) = NULL

test_that("Weighted network statistics match perfectly", { # 5 tests
  expect_equal(myNetStats[1], diameter(myG))
  expect_equal(myNetStats[2], mean(myDistances, na.rm = TRUE))
  expect_equal(myNetStats[2 + (1:N)], myBetweenness)
  expect_equal(myNetStats[2 + N + (1:N)], myCloseness)
  expect_equal(myNetStats[2 + 2*N + (1:N)], myEigenvector)
})

uTree = makeUnweighted(myTree)
uNetStats = computeNetworkStats(uTree, weight = TRUE, meanpath = TRUE, maxOnly = FALSE)
wNetStats = computeNetworkStats(myTree, weight = FALSE, meanpath = TRUE, maxOnly = FALSE)
uDF = as.data.frame(cbind(uTree$edge, weight = uTree$edge.length))
uG = igraph::graph_from_data_frame(uDF, directed = FALSE)
uN = igraph::gorder(uG)
uDistances = distances(uG)
diag(uDistances) = NA
uBetweenness = igraph::betweenness(uG, directed = FALSE)
unorder = order(as.numeric(names(uBetweenness)))
uBetweenness = uBetweenness[unorder]
uCloseness = igraph::closeness(uG)
uCloseness = uCloseness[unorder]
uEigenvector = igraph::eigen_centrality(uG, scale = FALSE)$vector
uEigenvector = uEigenvector[unorder]
names(uNetStats) = NULL
names(wNetStats) = NULL
names(uBetweenness) = NULL
names(uCloseness) = NULL
names(uEigenvector) = NULL

test_that("Unweighted network statistics match perfectly", { # 6 tests
  expect_identical(uNetStats, wNetStats)
  expect_equal(uNetStats[1], diameter(uG))
  expect_equal(uNetStats[2], mean(uDistances, na.rm = TRUE))
  expect_equal(uNetStats[2 + (1:N)], uBetweenness)
  expect_equal(uNetStats[2 + N + (1:N)], uCloseness)
  expect_equal(uNetStats[2 + 2*N + (1:N)], uEigenvector)
})
