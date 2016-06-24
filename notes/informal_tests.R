library(dpOMP)

set.seed(6251148)

G = 1
V = 1
K = 10
N = 10

d <- generate_data(X = diag(1), G=G, K=K, N=N)
out <- dpgmm(d$y, d$X, G, V, K, N, iter=10000)

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

####
G = 10000
V = 1
modelK = 50
trueK = 10
N = 9
d <- generate_data(X = diag(V), G=G, K=trueK, N=N)
d$z
out <- dpgmm(d$y, d$X, G, V, modelK, N, iter=10000)

#identification of true locations
dens = density(out$beta_g)
plot(dens)
abline(v=d$beta, lty=2)

# shrinkage
data_means = rowMeans(d$y)
post_means = rowMeans(out$beta_g[,,])
plot(data_means, post_means)
abline(0,1, col='red')

#correct pooling of information
hist(out$beta_g[,d$z==1,], prob=T, 30, main="correct pooling of information")
Gk = sum(d$z==1)
curve(dnorm(x, sum(d$xTy[d$z==1])/(Gk*N + 1), sd(d$y[d$z==2,])/sqrt(Gk*N + 1)), add=T, lty=2)
abline(v = mean(out$beta_g[,d$z==1,]))

#shrinkage?
opar = par(mfrow=c(4,4))
genes = sort(sample(G,16,replace=FALSE))
for (g in genes) {
hist(out$beta_g[,g,], prob=T, 30, main=paste("Gene=",g))
curve(dnorm(x, sum(d$xTy[g])/(N+1), sd(d$y[g,])/sqrt(N + 1)), add=T, lty=2)
abline(v = d$beta[d$z[g]])
}

#number of occupied clusters
occupied <- sapply(1:10000, function(i) sum(out$beta[1,,i] %in% out$beta_g[1,,i]))
hist(occupied, breaks = 1:50+.5)

#maximum index of significant pi
max_index <- sapply(1:10000, function(i) max(which(out$pi[,i]>.005)))
hist(max_index, breaks = 1:50+.5)
