#' Draw from full conditional for beta, sigma2 conditioning on Y=y, Z=z
#' @export
draw_theta <- function(z, yTy, xTy, xTx, G, K, V, N){
  .Call("draw_thetaR",
        as.integer(z),
        as.numeric(yTy),
        as.numeric(xTy),
        as.numeric(xTx),
        as.integer(G),
        as.integer(K),
        as.integer(V),
        as.integer(N),
        PACKAGE = "dpOMP")
}