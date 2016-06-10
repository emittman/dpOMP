#' generates a data set from the model
#' 
#' @export

generate_data <- function(X, n, G, K, beta=NULL){
  if(is.null(dim(X)) || dim(X)[1] != dim(X)[2]) stop("X must be a square matrix")
  V <- dim(X)[1]
  if(!is.null(beta)){
    if(length(beta) != V*G) stop("beta must be same dimension as X times K")
  } else{
    beta <- matrix(rnorm(K*V, 0, 2.5), V, K) 
  }
  Xexpand <- matrix(rep(X, each=n), V*n, V)
  groupMeans <- Xexpand %*% beta
  z <- sample(1:K, G, replace=T)
  y <- t(round(sapply(z, function(zz) rnorm(V*n, groupMeans[,zz])), 3))
  list(beta = beta,
       z = z,
       y=y,
       X = Xexpand,
       yTx = y%*%Xexpand,
       xTx = t(X)%*%X)
}