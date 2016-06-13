#' draw a simplex conditioned on current cluster allocations
#' 
#' @export
draw_pi <- function(Gk, K){
  .Call("draw_piR",
        as.numeric(Gk),
        as.integer(K),
        PACKAGE = "dpOMP")
}