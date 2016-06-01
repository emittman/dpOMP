#' Return first n square numbers
#' 
#' @param n how many d'ya want?
#'
#' @return it's pretty obvious
#'
#' @examples
#' n_squares(5)
#'
#' @export
n_squares <- function(n){
  .Call("n_squaresR",
        as.integer(n),
        PACKAGE = "dpOMP")
}