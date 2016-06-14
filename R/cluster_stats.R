#' Compute summary stats for groups
#' 
#' @param z group indicators
#' @param K max number of groups
#' 
#' @export
cluster_stats <- function(k, xTy, xTx, G, V, n, z){
  .Call("cluster_statsR",
        as.integer(k),
        as.numeric(xTy),
        as.numeric(xTx),
        as.integer(G),
        as.integer(V),
        as.integer(n),
        as.integer(z),
        PACKAGE = "dpOMP")
}
