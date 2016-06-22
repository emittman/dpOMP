context("Basic functionality")
# set.seed(61616)
# G <- 100
# V <- 2
# N <- 3
# X <- diag(V)
# 
# I <- 1000
# beta <- c(2, 5)
# d <- generate_data(X=X, N, G, K=1, beta=beta)
# 
# out <- with(d, dpgmm(data = y, design = diag(V), G, V, K=5, N, iter=I))
# 
# test_that("beta correct dimensions", dim(out$beta) == c(V, K, I))
# 
# test_that("pi correct dimensions", dim(out$pi) == c(K, I))
