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
    stats <- with(d, cluster_sums(k, yTy, xTy, G, V, z - 1))
    Gk <- sum(d$z == k+1)
    if(Gk > 1){
      yTyk <- with(d, sum(yTy[z == k+1]))
      xTyk <- with(d, rowSums(xTy[, z == k+1]))
    } else if(Gk ==1) {
      yTyk <- with(d, yTy[z == k+1])
      xTyk <- with(d, xTy[, z == k+1])
    } else {
      yTyk <- 0
      xTyk <- rep(0, V)
    }
    expect_equal(stats[[1]], Gk)
    expect_equal(stats[[2]], yTyk)
    expect_equal(stats[[3]], drop(xTyk))
#     stats[[1]]
#     Gk
#     stats[[2]]
#     yTyk
#     stats[[3]]
#     xTyk
  }
})