context("pi_prime")
V <- 3
x <- matrix(rnorm(V*V), V, V)
xTx <- crossprod(x)
y <- rnorm(V)
yTx <- y%*%x
pi1 <- runif(1,0,1)
pi2 <- runif(1,0,pi1)
beta1 <- rnorm(V)
beta2 <- rnorm(V)
sigma2 <- runif(1, 1, 2) * diag(V)

ccalc1 <- pi_prime(yTx, xTx, beta1, pi1, sigma2[1], V)
ccalc2 <- pi_prime(yTx, xTx, beta2, pi2, sigma2[1], V)

rcalc1 <- log(pi1) + -1/(2 * sigma2[1]) * (-2 * yTx %*% beta1 + t(beta1)%*%xTx%*%beta1)
rcalc2 <- log(pi2) + -1/(2 * sigma2[1]) * (-2 * yTx %*% beta2 + t(beta2)%*%xTx%*%beta2)

C <- ccalc1 - ccalc2
R <- rcalc1 - rcalc2
S <- log(pi1) - log(pi2) + mvtnorm::dmvnorm(y, x%*%beta1, sigma2, log = T) - mvtnorm::dmvnorm(y, x%*%beta2, sigma2, log=T)

test_that("check computation", {
  expect_equal(drop(C-R),  0)
})
test_that("matches 3rd party", {
  expect_equal(drop(R-S), 0)
})