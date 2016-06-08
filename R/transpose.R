#' transpose a matrix
#' @param matin M by N matrix or M*N vector
#' @param M integer
#' @param N integer
#' @export
transpose <- function(m, M, N){
  .Call("transposeR",
        as.numeric(m),
        as.integer(M),
        as.integer(N),
        PACKAGE = "dpOMP")
}