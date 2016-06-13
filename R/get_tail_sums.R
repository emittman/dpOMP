#' Get tail sums of a real-valued vector
#' 
#' @export

get_tail_sums <- function(vec){
  len <- length(vec)
  .Call("tail_sums",
        as.numeric(vec),
        as.integer(len),
        PACKAGE = "dpOMP")
} 