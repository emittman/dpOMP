context("cluster_stats")

set.seed(61016)
V <- 2
G <- 10
N <- 10
K <- 5
X <- matrix(c(1,1,1,-1), V, V)

d <- generate_data(X, N, G, K)

k <- as.numeric(names(which.max(table(d$z))))


stats <- with(d, cluster_stats(k, xTy, xTx, G, 2, N, z))

XTY <- with(d, rowSums(xTy[,z==k]))

Vinv <- d$xTx * stats[[1]] + diag(V)

V <- solve(Vinv)

bhat <- V %*% XTY

L <- chol(V)

test_that("estimate is correct", {
  expect_equal(stats[[2]], drop(bhat))
})

test_that("inverse is correct", {
  indic <- lower.tri(L, diag = T)
  expect_equal(indic*stats[[3]], L)
})

test_that("IG scale is correct", {
  expect_equal(stats[[4]], drop(t(bhat)%*%Vinv%*%bhat))
})