#' Return indicator that gene G is in cluster k
#' 
#' @export
get_k_indices <- function(z, k){
  G <- length(z)
  .Call("get_k_indicesR",
        as.integer(z),
        as.integer(G),
        as.integer(k),
        PACKAGE="dpOMP")
}