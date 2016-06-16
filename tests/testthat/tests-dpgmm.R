# context("Basic functionality")
set.seed(61616)
G <- 10
V <- 2
K <- 4
N <- 6
X <- diag(V)

d <- generate_data(X=X, N, G, K, pi=c(.9,.1/3,.1/3,.1/3))

out <- with(d, dpgmm(data = y, design = diag(V), G, V, K, N, iter=10))

data <- d$y
design <- X
Xexpand <- matrix(rep(design, each=N), V*N, V)
yTy <- sapply(1:G, function(g) data[g,]%*%data[g,])
yTy <- sum(yTy)
xTy <- t(Xexpand) %*% t(data)
xTx <- crossprod(design)*N
