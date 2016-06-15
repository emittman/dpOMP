#' Do MCMC for dpgmm model
#' 
#' @param data A column-major $\eqn{G \cdot VN}$ gene expressionvector/matrix
#' 
#' @param design A column-major $\eqn{V \cdot V}$ experimental design vector/matrix
#' 
#' @param G number of genes
#' 
#' @param K maximum number of clusters
#' 
#' @param V number of varieties/groups
#' 
#' @param N samples per variety/group
#' 
#' @export
dpgmm <- function(data, design, G, K, V, N){
  Xexpand <- matrix(rep(design, each=N), V*N, V)
  yTy <- sapply(1:G, function(g) data[g,]%*%data[g,])
  yTy <- sum(yTy)
  xTy <- t(Xexpand) %*% t(data)
  xTx <- crossprod(design)*N
  
  out <- .Call("dpgmmR",
               as.numeric(yTy),
               as.numeric(xTy),
               as.numeric(xTx),
               as.integer(G),
               as.integer(K),
               as.integer(V),
               as.integer(N),
               PACKAGE = "dpOMP")
  names(out) <- c("beta", "pi")
  return(out)
}