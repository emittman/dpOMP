rstar <- function(n, p1, p2, lambda1, lambda2){
	A <- 2*pi*runif(n,0,1)
	upper <- p1*sin(lambda1*A) + p2*cos(lambda2*A)
	r <- sign(upper) * runif(n, 0, abs(upper))
	out <- data.frame(xout = r*cos(A), yout = r*sin(A))
	return(out)
}

#V <- 2 
G <- 500
N <- 3
sigma2 <- .1

beta_sim = rstar(G, 1, 1, 2, 2)
plot(beta_sim)

beta <- apply(beta_sim, 1, identity)
X <- diag(2) %x% rep(1, N)
y <- matrix(rnorm(G*2*N, X%*%beta, sigma2), G, 2*N, byrow=T)

out <- dpgmm(data=y, design=X, G, 2, 100, N, 10000)
hist(out$sigma2, 30)

par(mfrow=c(4, 4))
for(g in 1:16){
  hist(out$beta_g[1,g,], 30)
  abline(v=beta[1,g])
}

library(ggplot2)
b1 <- as.numeric(out$beta_g[1,,])
b2 <- as.numeric(out$beta_g[2,,])
b <- data.frame(xout=b1, yout=b2)
ggplot(b, aes(xout, yout)) + geom_hex(bins=20) + geom_point(data=beta_sim) +
  scale_fill_continuous(trans="log")
