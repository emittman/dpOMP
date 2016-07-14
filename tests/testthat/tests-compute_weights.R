context("compute_weights")
require(plyr)
num_tests <- 10
set.seed(621621)
data_sets <- llply(1:num_tests, function(i){
  G = rpois(1, 9.5)+1
  K = rpois(1, 2.5)+1
  V = rpois(1, 1.5)+1
  x = matrix(rnorm(V*V),V,V)
  d <- generate_data(X=x, N=3, G=G, K=K)
  return(d)
})


  outputs <- llply(1:num_tests, function(i){
    data <- data_sets[[i]]
    with(data, compute_weights(yTy, xTy, xTx, beta, pi, sigma2, G, K, V))
  })

test_that("output size correct", {
  for(i in 1:num_tests){
    #print(length(outputs[[i]]) == with(data_sets[[i]], G*K))
    expect_equal(length(outputs[[i]]), with(data_sets[[i]], G * K))
  }
})

pi_primeR <- function(yTy, xTy, xTx, beta, pi, sigma2, V){
  log(pi) -log(sigma2)  + - 1.0/(2.0*sigma2) * (yTy + t(beta)%*% xTx %*% beta - 2 * xTy %*% beta)
}

compute_weightsR <- function(yTy, xTy, xTx, beta, pi, sigma2, G, K, V){
  grid <- expand.grid(k=1:K, g=1:G)
  daply(grid, .(k,g), function(cell) {
    with(cell, pi_primeR(yTy[g], xTy[(V*(g-1)+1):(V*g)], xTx, beta[(V*(k-1)+1):(V*k)], pi[k], sigma2[k], V=V))
         }
  )}

outputsR <- llply(1:num_tests, function(i){
  data <- data_sets[[i]]
  out <- with(data, compute_weightsR(yTy, xTy, xTx, beta, pi, sigma2, G, K, V))
  dim(out) <- NULL
  out
  })

test_that("output is correct", {
  for(i in 1:num_tests){
    #print(all.equal(outputs[[i]], outputsR[[i]]))
    expect_equal(outputs[[i]], outputsR[[i]])
  }
})

##?delete
# #   data2 <- list(
# #     xTy = c(1,1),
# #     xTx = diag(2),
# #     beta = c(1,-.2),
# #     pi = 1,
# #     G = 1,
# #     K = 1,
# #     n = 1,
# #     V = 2
# #   )
# #   out2 <- with(data1, compute_weights(xTy, xTx, beta, pi, G, K, n, V))
# #   pp <- with(data1, pi_prime(xTy, xTx, beta, pi, 1.0, V))
#   pi_primeR <- function(xTy, xTx, beta, pi, sigma2, V){
#     log(pi)  + - 1.0/(2.0*sigma2) * (t(beta)%*% xTx %*% beta - 2 * xTy %*% beta)
#   }
#   compute_weightsR <- function(xTy, xTx, beta, pi, sigma2, G, K, n, V){
#     g <- expand.grid(g=1:G, k=1:K)
#   }
# #   with(data2, pi_primeR(xTy, xTx, beta, pi, 1.0, V))
# #   with(data2, pi_prime(xTy, xTx, beta, pi, 1.0, V))

