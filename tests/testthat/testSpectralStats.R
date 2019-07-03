library(treeCentrality)
library(ape)
library(igraph)
library(RPANDA)
context("Spectral statistics on a tree")

myTree = createTestTree()
mySpecStats = computeSpectralStats(myTree, weight = TRUE, adj = c(FALSE, TRUE),
    norm = c(FALSE, TRUE), dist = c(FALSE, TRUE), full = c(FALSE, TRUE), maxOnly = FALSE, unitMean = FALSE)
myDF = as.data.frame(cbind(myTree$edge, weight = myTree$edge.length))
myG = igraph::graph_from_data_frame(myDF, directed = FALSE)
Spectrum = eigen(as_adjacency_matrix(myG, attr = "weight", sparse = FALSE))$values
LapSpectrum = eigen(laplacian_matrix(myG, normalized = FALSE, sparse = FALSE))$values
NormLapSpectrum = eigen(laplacian_matrix(myG, normalized = TRUE, sparse = FALSE))$values
spectR1 = RPANDA::spectR(myTree, "standard")$eigenvalues
# spectR2 = spectR(myTree, "normal")$eigenvalues - this function uses asymmetric matrices

test_that("Weighted spectral statistics match perfectly", { # 4 tests
  expect_equal(mySpecStats[[1]], LapSpectrum)
  # expect_equal(mySpecStats[[2]], DistLapSpectrum) - skipping this test for now
  expect_equal(mySpecStats[[3]], c(spectR1, 0))
  expect_equal(mySpecStats[[4]], NormLapSpectrum)
  # expect_equal(mySpecStats[[5]], NormDistLapSpectrum) - skipping this test for now
  # expect_equal(mySpecStats[[6]], spectR2) - skipping this test for now
  expect_equal(mySpecStats[[7]], Spectrum)
  # expect_equal(mySpecStats[[8]], DistSpectrum) - skipping this test for now
  # expect_equal(mySpecStats[[9]], FullDistSpectrum) - skipping this test for now
})

uTree = makeUnweighted(myTree)
uSpecStats = computeSpectralStats(uTree, weight = TRUE, adj = c(FALSE, TRUE),
    norm = c(FALSE, TRUE), dist = c(FALSE, TRUE), full = c(FALSE, TRUE), maxOnly = FALSE, unitMean = FALSE)
names(uSpecStats) = NULL
wSpecStats = computeSpectralStats(myTree, weight = FALSE, adj = c(FALSE, TRUE),
    norm = c(FALSE, TRUE), dist = c(FALSE, TRUE), full = c(FALSE, TRUE), maxOnly = FALSE)
names(wSpecStats) = NULL
uDF = as.data.frame(cbind(uTree$edge, weight = uTree$edge.length))
uG = igraph::graph_from_data_frame(uDF, directed = FALSE)
uSpectrum = eigen(as_adjacency_matrix(uG, attr = "weight", sparse = FALSE))$values
uLapSpectrum = eigen(laplacian_matrix(uG, normalized = FALSE, sparse = FALSE))$values
uNormLapSpectrum = eigen(laplacian_matrix(uG, normalized = TRUE, sparse = FALSE))$values
uspectR1 = RPANDA::spectR(uTree, "standard")$eigenvalues
# uspectR2 = spectR(uTree, "normal")$eigenvalues  - this function uses asymmetric matrices

test_that("Unweighted spectral statistics match perfectly", { # 5 tests
  expect_identical(uSpecStats, wSpecStats)
  expect_equal(uSpecStats[[1]], uLapSpectrum)
  # expect_equal(uSpecStats[[2]], uDistLapSpectrum) - skipping this test for now
  expect_equal(uSpecStats[[3]], c(uspectR1, 0))
  expect_equal(uSpecStats[[4]], uNormLapSpectrum)
  # expect_equal(uSpecStats[[5]], uNormDistLapSpectrum) - skipping this test for now
  # expect_equal(uSpecStats[[6]], uspectR2) - skipping this test for now
  expect_equal(uSpecStats[[7]], uSpectrum)
  # expect_equal(uSpecStats[[8]], uDistSpectrum) - skipping this test for now
  # expect_equal(uSpecStats[[9]], uFullDistSpectrum) - skipping this test for now
})
