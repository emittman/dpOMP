#' Get component weights for all genes
#' 
#' @param yTx numeric vector of length G*V
#' @param xTx numeric vector of length V*V
#' @param beta numeric vector of length K*V
#' @param pi numeric vector of length K
#' @param sigma2 double
#' @param G int
#' @param K int
#' @param n int
#' @param V int
#' @export
compute_weights <- function(yTy, yTx, xTx, beta, pi, sigma2, G, K, V){
  out <- .Call("compute_weightsR",
               as.numeric(yTy),
               as.numeric(yTx),
               as.numeric(xTx),
               as.numeric(beta),
               as.numeric(pi),
               as.numeric(sigma2),
               as.integer(G),
               as.integer(K),
               as.integer(V),
               PACKAGE = "dpOMP")
  return(out)
}
