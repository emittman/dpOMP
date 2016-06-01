context("weights")

test_that("weights returned", {
  data1 <- list(
    yTx  = rnorm(10*2),
    xTx  = diag(2),
    beta = rnorm(3*2),
    pi   = rep(1/3, 3),
    G    = 10,
    K    = 3,
    n    = 1,
    V    = 2
  )
  out1 <- with(data1, compute_weights(yTx, xTx, beta, pi, G, K, n, V))
  expect_equal(length(out1), data1$G*data1$K)
  data2 <- list(
    yTx = c(1,1),
    xTx = diag(2),
    beta = c(1,-.2),
    pi = 1,
    G = 1,
    K = 1,
    n = 1,
    V = 2
  )
  out2 <- with(data1, compute_weights(yTx, xTx, beta, pi, G, K, n, V))
  pp <- with(data1, pi_prime(yTx, xTx, beta, pi, 1.0, V))
  pi_primeR <- function(yTx, xTx, beta, pi, sigma2, V){
    log(pi)  + - 1.0/(2.0*sigma2) * (t(beta)%*% xTx %*% beta - 2 * yTx %*% beta)
  }
  compute_weightsR <- function(yTx, xTx, beta, pi, sigma2, G, K, n, V){
    g <- expand.grid(g=1:G, k=1:K)
  }
  with(data2, pi_primeR(yTx, xTx, beta, pi, 1.0, V))
  with(data2, pi_prime(yTx, xTx, beta, pi, 1.0, V))
})
