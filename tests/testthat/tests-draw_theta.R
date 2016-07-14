context("draw theta")
set.seed(61516)
V <- 2
X <- diag(V)
G <- 6
K <- 4
N <- 10

d <- generate_data(X=X, N=N, G=G, K=K, beta=matrix(rep(c(1,-1), K), V, K), pi=rep(1/K, K))
with(d, draw_theta(z-1, yTy, xTy, xTx, G, K, V, N))
d$z-1
table(d$z-1)
