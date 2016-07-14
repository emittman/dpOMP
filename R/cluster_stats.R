#' Compute sums for a group
#' 
#' @param k which group (zero indexed)
#' @param yTy numeric vector length G
#' @param xTy numeric matrix V by G
#' @param G integer
#' @param V integer
#' @param z group indicators
#' 
#' @export
cluster_sums <- function(k, yTy, xTy, G, V, z){
  .Call("cluster_sumsR",
        as.integer(k),
        as.numeric(yTy),
        as.numeric(xTy),
        as.integer(G),
        as.integer(V),
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
#' @param yTyk squared norm of data associated with cluster k
#' @param xTyk crossproduct of data associated with cluster k and design matrix
#' @param xTx  crossproduct of design matrix
#' @param betak current value of beta parameter for cluster k
#' @param Gk number of genes associated with cluster k
#' @export
calculate_IGscale <- function(yTyk, xTyk, xTx, betak, Gk){
  out <- .Call("calculate_IGscaleR",
               as.numeric(yTyk),
               as.numeric(xTyk),
               as.numeric(xTx),
               as.numeric(betak),
               as.numeric(Gk),
               as.integer(length(xTyk)))
  return(out)
}