#' dot product of two vectors
#' 
#' @export
dot_prod <- function(x, y){
  .Call("ddotR",
        as.integer(length(x)),
        as.numeric(x),
        as.numeric(y),
        PACKAGE = "dpOMP")
}