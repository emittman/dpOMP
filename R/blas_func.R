#' BLAS implemented daxpy routine
#' 
#' 
linear_comb_vec <- function(alpha, x, y){
  .Call("linear_comb_vecR",
        as.integer(length(x)),
        as.numeric(alpha),
        as.numeric(x),
        as.numeric(y),
        PACKAGE = "dpOMP"
        )
}

#' BLAS implemented dsymv routine
#' 
#' 
dsymv <- function(alpha, A, x, beta, y){
  .Call("dsymvR",
        as.integer(length(x)),
        as.numeric(alpha),
        as.numeric(A),
        as.numeric(x),
        as.numeric(beta),
        as.numeric(y),
        PACKAGE = "dpOMP"
  )
}

#' Computes xTAx using BLAS routines dsymv and ddot (through wrappers)
#' 
#' @export
quad_form <- function(x, A){
  .Call("quad_formR",
        as.integer(length(x)),
        as.numeric(x),
        as.numeric(A),
        PACKAGE = "dpOMP")
}
