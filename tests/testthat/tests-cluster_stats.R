context("cluster_stats")
G = 20
V = 3
n = 4
X = matrix(rep(c(1, 1, 1,
                 1, -1, 1,
                 -1, 1, 1), each=n), 3*n, 3)
xTx <- t(X)%*%X



