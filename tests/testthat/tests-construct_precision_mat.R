context("construct precision matrix")
V <- 5
Gk <- 3
X <- matrix(rnorm(V*V), V, V)
xTx <- crossprod(X)
sigma2 <- 1/rgamma(1, 1, 1)
lambda2 <- 1/rgamma(1, 1, 1)
SinvR <- xTx*Gk + diag(rep(sigma2/lambda2, V))
SinvC <- construct_precision_mat(xTx, Gk, sigma2, lambda2, V)
test_that("Prec. matrix is correct",{
  expect_equal(as.numeric(SinvR), as.numeric(SinvC))
})