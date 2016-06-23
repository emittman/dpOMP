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
dpgmm <- function(data, design, G, V, K, N, iter){
  yTy <- sapply(1:G, function(g) data[g,]%*%data[g,])
  yTy <- sum(yTy)
  xTy <- t(design) %*% t(data)
  xTx <- crossprod(design)
  
  out <- .Call("dpgmmR",
               as.numeric(yTy),
               as.numeric(xTy),
               as.numeric(xTx),
               as.integer(G),
               as.integer(V),
               as.integer(K),
               as.integer(N),
               as.integer(iter),
               PACKAGE = "dpOMP")
  names(out) <- c("beta", "pi", "beta_g", "sigma2")
  out$beta <- array(out$beta, dim=c(V, K, iter))
  out$pi <- array(out$pi, dim=c(K, iter))
  out$beta_g <- array(out$beta_g, dim=c(V, G, iter))
  return(out)
}