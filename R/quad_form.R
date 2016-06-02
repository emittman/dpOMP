# #' Compute xTAx
# #' 
# #' @param A real-valued symmetric matrix
# #' @param x real-valued vector with dimension compatible with A
# #' @param dim integer giving dimensions of square matrix A
# #' @return numeric vector of length 1
# #' @export
# quad_form <- function(A, x, dim){
#   out <- .Call("quad_formR",
#                as.numeric(A),
#                as.numeric(x),
#                as.integer(dim),
#                PACKAGE = "dpOMP")
#   return(out)
# }
