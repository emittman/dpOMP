require(plyr)
require(MCMCpack)
context("weights")
num_tests <- 10
set.seed(621621)
data_sets <- llply(1:num_tests, function(i){
  G = rpois(1, 9.5)+1
  K = rpois(1, 2.5)+1
  V = rpois(1, 1.5)+1
  x = matrix(rnorm(V*V),V,V)
  list(
    yTx  = rnorm(G*V),
    xTx  = t(x)%*%x,
    beta = rnorm(K*V),
    pi   = rdirichlet(1, rep(1,K)),
    G    = G,
    K    = K,
    n    = 1,
    V    = V
  )
})


  outputs <- llply(1:num_tests, function(i){
    data <- data_sets[[i]]
    with(data, compute_weights(yTx, xTx, beta, pi, G, K, n, V))
  })

test_that("output size correct", {
  for(i in 1:num_tests)
    expect_equal(length(outputs[[i]]), with(data_sets[[i]], G*K))
})

pi_primeR <- function(yTx, xTx, beta, pi, sigma2=1, V){
  log(pi)  + - 1.0/(2.0*sigma2) * (t(beta)%*% xTx %*% beta - 2 * yTx %*% beta)
}
compute_weightsR <- function(yTx, xTx, beta, pi, G, K, n, V){
  grid <- expand.grid(k=1:K, g=1:G)
  daply(grid, .(k,g), function(cell) {
    with(cell, pi_primeR(yTx[(V*(g-1)+1):(V*g)], xTx, beta[(V*(k-1)+1):(V*k)], pi[k], V=V))
         }
  )}

outputsR <- llply(1:num_tests, function(i){
  data <- data_sets[[i]]
  out <- with(data, compute_weightsR(yTx, xTx, beta, pi, G, K, n, V))
  dim(out) <- NULL
  out
  })

test_that("output is correct", {
  for(i in 1:num_tests)
    expect_equal(outputs[[i]], outputsR[[i]])
})
#   data2 <- list(
#     yTx = c(1,1),
#     xTx = diag(2),
#     beta = c(1,-.2),
#     pi = 1,
#     G = 1,
#     K = 1,
#     n = 1,
#     V = 2
#   )
#   out2 <- with(data1, compute_weights(yTx, xTx, beta, pi, G, K, n, V))
#   pp <- with(data1, pi_prime(yTx, xTx, beta, pi, 1.0, V))
  pi_primeR <- function(yTx, xTx, beta, pi, sigma2, V){
    log(pi)  + - 1.0/(2.0*sigma2) * (t(beta)%*% xTx %*% beta - 2 * yTx %*% beta)
  }
  compute_weightsR <- function(yTx, xTx, beta, pi, sigma2, G, K, n, V){
    g <- expand.grid(g=1:G, k=1:K)
  }
#   with(data2, pi_primeR(yTx, xTx, beta, pi, 1.0, V))
#   with(data2, pi_prime(yTx, xTx, beta, pi, 1.0, V))

