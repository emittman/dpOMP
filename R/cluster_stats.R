#' Compute sums for groups
#' 
#' @param z group indicators
#' @param K max number of groups
#' 
#' @export
cluster_sums <- function(k, xTy, G, V, N, z){
  .Call("cluster_sumsR",
        as.integer(k),
        as.numeric(xTy),
        as.integer(G),
        as.integer(V),
        as.integer(N),
        as.integer(z),
        PACKAGE = "dpOMP")
}

#' Construct precision matrix for full conditional distribution for beta_k
#' 
#' @param xTx Crossproduct of design matrix
#' @param Gk Number of genes associated with the cluster
#' @param sigma2 Data variance
#' @param lambda2 Prior variance
#' @param V dimension of beta
#' @export
construct_precision_mat <- function(xTx, Gk, sigma2, lambda2, V){
  out <- .Call("construct_precision_matR",
               as.numeric(xTx),
               as.numeric(Gk),
               sigma2,
               lambda2,
               as.integer(V))
  return(as.matrix(out, V, V))
}

#' Calculate cluster contribution to IG scale paramter of full conditional for sigma2
#' 
#' @param xTyk Crossproduct of data associated with cluster k and design matrix
#' @param xTx Crossproduct of design matrix
#' @param betak Current value of beta parameter for cluster k
#' @param Gk Number of genes associated with cluster k
#' @export
increment_IGscale <- function(xTyk, xTx, betak, Gk){
  out <- .Call("increment_IGscaleR",
               as.numeric(xTyk),
               as.numeric(xTx),
               as.numeric(betak),
               as.numeric(Gk),
               as.integer(length(xTyk)))
  return(out)
}