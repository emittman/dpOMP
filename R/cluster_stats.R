#' Compute summary stats for groups
#' 
#' @param z group indicators
#' @param K max number of groups
#' 
#' @export
cluster_stats <- function(z, yTx, K, V){
  .Call("cluster_statsR",
        as.integer(z),
        as.numeric(yTx),
        as.integer(length(z)),
        as.integer(K),
        as.integer(V),
        PACKAGE = "dpOMP")
}