#' Compute sums for groups
#' 
#' @param z group indicators
#' @param K max number of groups
#' 
#' @export
cluster_sums <- function(k, xTy, G, V, N, z){
  .Call("cluster_sumsR",
        as.integer(k),
        as.numeric(xTy),
        as.integer(G),
        as.integer(V),
        as.integer(N),
        as.integer(z),
        PACKAGE = "dpOMP")
}
