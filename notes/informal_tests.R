G = 1
V = 1
K = 10
N = 10

d <- generate_data(X = diag(1), G=G, K=K, N=N)

w <- compute_weights(d$yTx, d$xTx, d$beta, d$pi, 1.0, G, K, V)
draw_z(weights = w, G=G, K=K)

out <- dpgmm(d$y, d$X, G, V, K, N, iter=1)
