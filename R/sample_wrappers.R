#' Draw random inverse gamma variates
#' @export
rinvgamma <- function(n, a, b){
  .Call("rinvgammaR",
        as.integer(n),
        a, b,
        PACKAGE = "dpOMP")
}