extern "C" SEXP send_chain(SEXP yTy, SEXP xTy, SEXP xTx, SEXP lambda2R, SEXP alphaR, SEXP aR, SEXP bR, SEXP G, SEXP K, SEXP V, SEXP N){

  set.seed(61516)
  V <- 2
  X <- diag(V)
  G <- 6
  K <- 4
  N <- 10
  
  d <- generate_data(X=X, N=N, G=G, K=K, pi=rep(1/K, K), sigma2=rep(1, K))
  
  with(d, .Call("send_chain",
                yTy,
                xTy,
                xTx,
                1, 1, 1, 1,
                as.integer(G),
                as.integer(K),
                as.integer(V),
                as.integer(N)))
  
  w <- with(d, compute_weights(yTy, xTy, xTx, beta=beta, pi, sigma2, G, K, V))
  dim(w) <- c(K, G)
  i=2
  d$z[i]
  d$beta[,d$z[i]]
  d$xTy[,i]/(d$N)
  d$beta
  w[,i]
  