G = 1
V = 1
K = 10
N = 10

d <- generate_data(X = diag(1), G=G, K=K, N=N)
out <- dpgmm(d$y, d$X, G, V, K, N, iter=1000)
#hist out$beta_g should match N(d$xTy/11, sd(d$y)/sqrt(11))
#hist out$sigma2 should match IG(5, d$xTy^2/11)


#compare true weights to proportion sampled by stick-breaking; some bias toward lower indices due to prior
w <- compute_weights(d$yTx, d$xTx, d$beta, d$pi, 1.0, G, K, V)
w <- log(exp(w)/sum(exp(w)))
s <- sapply(1:1e5, function(i) draw_z(weights = w, G=G, K=K))
prop_tbl <- rbind(table(s)/1e5,
      round(exp(w)[as.integer(names(table(s))) + 1], 5),
      round(table(s)/1e5- exp(w)[as.integer(names(table(s))) + 1],5))
row.names(prop_tbl) <- c("actual","expected","difference")
prop_tbl
