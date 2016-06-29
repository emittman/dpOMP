rstar <- function(n, p1, p2, lambda1, lambda2){
	A <- 2*pi*runif(n,0,1)
	upper <- p1*sin(lambda1*A) + p2*cos(lambda2*A)
	r <- sign(upper) * runif(n, 0, abs(upper))
	out <- data.frame(xout = r*cos(A), yout = r*sin(A))
	return(out)
}

#???
plot(rstar(1000, 2, 2, 5, 4))

plot(rstar(1e4, 1, 1, 3, 5))

beta_sim = rstar(10000, 0, 1, 3, 3)
plot(beta_sim)

