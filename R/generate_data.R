#' generates a data set from the model
#' 
#' @param X V by V numeric design matrix
#' @param N number of replicates per variety
#' @param G number of genes
#' @param K number of clusters
#' @param beta V by K matrix of location parameters
#' @param pi K-vector of cluster weights
#' 
#' @export

generate_data <- function(X, N, G, K, beta=NULL, pi=NULL, sigma2=1.0){
  if(is.null(dim(X)) || dim(X)[1] != dim(X)[2]) stop("X must be a square matrix")
  V <- dim(X)[1]
  if(!is.null(beta)){
    if(length(beta) != V*K) stop("beta must have length K times rank(X)")
  } else{
    beta <- matrix(rnorm(K*V, 0, 2.5), V, K) 
  }
  if(!is.null(pi)){
    if(length(pi) != K) stop("pi must have dimension K")
  } else{
    pi <- MCMCpack::rdirichlet(1, rep(1, K)) 
  }
  Xexpand <- matrix(rep(X, each=N), V*N, V)
  groupMeans <- Xexpand %*% beta
  z <- sample(1:K, G, replace=T, prob = pi)
  y <- t(round(sapply(z, function(zz) rnorm(V*N, groupMeans[,zz], sigma2)), 3))
  list(G = G,
       V = V,
       N = N,
       K = K,
       beta = beta,
       pi = pi,
       z = z,
       y=y,
       X = Xexpand,
       xTy = t(Xexpand)%*%t(y),
       xTx = t(Xexpand)%*%Xexpand)
}