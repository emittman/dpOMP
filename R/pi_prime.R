#' Update component weight for given data
#' 
#' @param yTx numeric vector of length V
#' @param xTx numeric vector of length V*V
#' @param beta numeric vector of length V
#' @param pi positive real number
#' @param sigma2 positive real number
#' @param V int
#' @export
pi_prime <- function(yTx, xTx, beta, pi, sigma2, V){
  if(pi<0 | pi>1) error("pi must be between 0 and 1")
  if(sigma2<= 0) error("sigma2 must be positive")
  out <- .Call("pi_primeR",
               as.numeric(yTx),
               as.numeric(xTx),
               as.numeric(beta),
               as.numeric(pi),
               as.numeric(sigma2),
               as.integer(V),
               PACKAGE="dpOMP")
  return(out)
}

