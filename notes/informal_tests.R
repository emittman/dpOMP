G = 1
V = 1
K = 10
N = 10

d <- generate_data(X = diag(1), G=G, K=K, N=N)
out <- dpgmm(d$y, d$X, G, V, K, N, iter=1000)
#hist out$beta_g should match N(d$xTy/11, sd(d$y)/sqrt(11))
hist(out$beta_g, prob=T, 30)
curve(dnorm(x, d$xTy/11, sd(d$y)/sqrt(11)), add=T, lty=2)
#hist out$sigma2 should match IG(N/2, (sum(d$y^2) - d$xTy^2/11)/2)
hist(out$sigma2, prob=T, 50)
curve(MCMCpack::dinvgamma(x, N/2, (sum(d$y^2) - d$xTy^2/11)/2), add=T, lty=2)

#compare true weights to proportion sampled by stick-breaking; some bias toward lower indices due to prior
w <- compute_weights(d$yTx, d$xTx, d$beta, d$pi, 1.0, G, K, V)
w <- log(exp(w)/sum(exp(w)))
s <- sapply(1:1e5, function(i) draw_z(weights = w, G=G, K=K))
prop_tbl <- rbind(table(s)/1e5,
      round(exp(w)[as.integer(names(table(s))) + 1], 5),
      round(table(s)/1e5- exp(w)[as.integer(names(table(s))) + 1],5))
row.names(prop_tbl) <- c("actual","expected","difference")
prop_tbl
