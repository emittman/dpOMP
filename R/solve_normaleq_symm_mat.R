#' Solves Ax = b, replacing input b with solution and A with cholesky factor
#' 
#' @export
solve_normaleq_symm_mat <- function(n, A, B){
  .Call("solve_normaleq_symm_matR",
        as.integer(n),
        as.numeric(A),
        as.numeric(B),
        PACKAGE = "dpOMP")
}