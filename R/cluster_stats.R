#' Compute summary stats for groups
#' 
#' @param z group indicators
#' @param K max number of groups
#' 
#' @export
cluster_stats <- function(z, K){
  .Call("cluster_statsR",
        as.integer(z),
        as.integer(length(z)),
        as.integer(K),
        PACKAGE = "dpOMP")
}