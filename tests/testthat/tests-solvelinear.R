context("solve linear equations")

set.seed(61016)
V <- 2
X <- matrix(c(1,1,1,-1),V,V)
G <- 10
n <- 3
K <- 5
d <- generate_data(X, n, G, K)

k <- as.numeric(names(which.max(table(d$z))))

stats <- with(d, cluster_stats(k, yTx, xTx, G, 2, n, z))

YTX <- with(d, colSums(yTx[z==k,]))

V <- solve(d$xTx * stats[[1]] * n + diag(V))

bhat <- V %*% YTX

L <- chol(V)

test_that("estimate is correct", {
  expect_equal(stats[[2]], drop(bhat))
})

test_that("inverse is correct", {
  indic <- lower.tri(L, diag = T)
  expect_equal(indic*stats[[3]], L)
})