#' Get a multinomial draw given specified weights
#' 
#' @param weights
#' 
#' @export
rcategorical <- function(weights){
  .Call("rcategoricalR",
        as.numeric(weights),
        as.integer(length(weights)),
        PACKAGE = "dpOMP")
}


#' Subtract max element and exponentiate lp vector
#' @export
norm_exp_lp <- function(weights, seed){
  .Call("norm_exp_lpR",
        as.numeric(weights),
        as.integer(length(weights)),
        PACKAGE = "dpOMP")
}

