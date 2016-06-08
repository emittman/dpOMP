#' Solves Ax = b, replacing input b with solution and A with cholesky factor
#' 
#' @export
dposv <- function(n, A, B){
  .Call("dposvR",
        as.integer(n),
        as.numeric(A),
        as.numeric(B),
        PACKAGE = "dpOMP")
}