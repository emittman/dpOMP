#' generates a data set from the model
#' 
#' @param X V by V numeric design matrix
#' @param N number of replicates per variety
#' @param G number of genes
#' @param K number of clusters
#' @param beta V by K matrix of cluster locations
#' @param pi K-vector of cluster weights
#' @param sigma2 K-vector of cluster variances
#' 
#' @export

generate_data <- function(X, N, G, K, beta=NULL, pi=NULL, sigma2=NULL){
  if(is.null(dim(X)) || dim(X)[1] != dim(X)[2]) stop("X must be a square matrix")
  V <- dim(X)[1]
  if(!is.null(beta)){
    if(length(beta) != V*K) stop("beta must have length K times rank(X)")
  } else{
    #default lambda2 = 2.5^2
    beta <- matrix(rnorm(K*V, 0, 2.5), V, K) 
  }
  if(!is.null(pi)){
    if(length(pi) != K) stop("pi must have dimension K")
  } else{
    pi <- MCMCpack::rdirichlet(1, rep(1, K)) 
  }
  if(!is.null(sigma2)){
    if(length(sigma2) !=K) stop("sigma2 must have dimension K")
  } else{
    #default a = b = 1
    sigma2 <- MCMCpack::rinvgamma(K, 1, 1)
  }
  Xexpand <- matrix(rep(X, each=N), V*N, V)
  groupMeans <- Xexpand %*% beta
  z <- sample(1:K, G, replace=T, prob = pi)
  y <- t(round(sapply(z, function(zz) rnorm(V*N, groupMeans[,zz], sqrt(sigma2[zz]))), 3))
  list(G = G,
       V = V,
       N = N,
       K = K,
       beta = beta,
       pi = pi,
       sigma2 = sigma2,
       z = z,
       y=y,
       X = Xexpand,
       yTy = sapply(1:G, function(g) crossprod(y[g,])),
       xTy = t(Xexpand)%*%t(y),
       xTx = t(Xexpand)%*%Xexpand)
}