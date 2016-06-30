context("cluster_sums")

set.seed(61016)
V <- 2
G <- 10
N <- 10
K <- 5
X <- matrix(c(1,1,1,-1), V, V)

d <- generate_data(X, N, G, K)

# k <- as.numeric(names(which.max(table(d$z))))

test_that("cluster sums are correct", {

  for(k in 0:(K-1)){
    stats <- with(d, cluster_sums(k, xTy, G, V, N, z))
    Gk <- sum(d$z == k+1)
    xTyk <- with(d, rowSums(xTy[, z == k+1]))
    expect_equal(stats[[1]], Gk)
    expect_equal(stats[[2]], xTyk)
  }
})



# Vinv <- d$xTx * stats[[1]] + diag(V)
# 
# V <- solve(Vinv)
# 
# bhat <- V %*% XTY
# 
# L <- chol(V)

# test_that("estimate is correct", {
#   expect_equal(stats[[2]], drop(bhat))
# })
# 
# test_that("inverse is correct", {
#   indic <- lower.tri(L, diag = T)
#   expect_equal(indic*stats[[3]], L)
# })
# 
# test_that("IG scale is correct", {
#   expect_equal(stats[[4]], drop(t(bhat)%*%Vinv%*%bhat))
# })