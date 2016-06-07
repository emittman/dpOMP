#' Draw z for all genes given the computed weights
#' 
#' @param weights row-major vector of unnormalized weights for G genes
#' @param G number of genes
#' @param K number of categories
#' 
#' @export
draw_z <- function(weights, G , K){
  .Call("draw_zR",
        as.numeric(weights),
        as.integer(G),
        as.integer(K),
        PACKAGE = "dpOMP")
}