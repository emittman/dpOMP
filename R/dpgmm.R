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
  out <- .Call("dpgmmR",
               as.numeric(data),
               as.numeric(design),
               as.integer(G),
               as.integer(K),
               as.integer(V),
               as.integer(N),
               PACKAGE = "dpOMP")
  names(out) <- c("beta", "pi")
  return(out)
}