#' Run MCMC for DPGMMM model, utilizing pilot chain
#' 
#' @param data A column-major $\eqn{G \cdot VN}$ gene expressionvector/matrix
#' @param design A column-major $\eqn{V \cdot V}$ experimental design vector/matrix
#' @param G number of genes
#' @param K maximum number of clusters
#' @param V number of varieties/groups
#' @param N samples per variety/group
#' @param iter number of iterations for main chain
#' @param init_iter number of iterations for pilot chain
#' 
#' @export

dpgmm_init <- function(data, design, G, V, K, N, iter, init_iter){
  # run an initial chain
  yTy <- sapply(1:G, function(g) data[g,]%*%data[g,])
  yTy <- sum(yTy)
  xTy <- t(design) %*% t(data)
  xTx <- crossprod(design)
  
  init_out <- .Call("dpgmmR",
               as.numeric(yTy),
               as.numeric(xTy),
               as.numeric(xTx),
               as.integer(G),
               as.integer(V),
               as.integer(K),
               as.integer(N),
               as.integer(init_iter),
               PACKAGE = "dpOMP")
  names(init_out) <- c("beta", "pi", "beta_g", "sigma2")
  init_out$beta <- array(out$beta, dim=c(V, K, iter))
  init_out$pi <- array(init_out$pi, dim=c(K, iter))
  init_out$beta_g <- array(init_out$beta_g, dim=c(V, G, iter))
  init <- list()
  # move most relevant atoms to front of vector
  indices <- which(init_out$beta[1,,init_iter] %in% init_out$beta_g[1,,init_iter])
  anti_indices <- (1:K)[!(1:K %in% indices)]
  indices <- indices[order(-init_out$pi[indices,init_iter])]
  init$beta <- init_out$beta_g[,c(indices, anti_indices),init_iter]
  init$pi <- init_out$pi[c(indices, anti_indices),init_iter]
  init$sigma2 <- init_out$sigma2[init_iter]
  
  out <- .Call("dpgmm_initR",
               as.numeric(yTy),
               as.numeric(xTy),
               as.numeric(xTx),
               as.numeric(init$beta),
               as.numeric(init$pi),
               as.numeric(init$sigma2),
               as.integer(G),
               as.integer(V),
               as.integer(K),
               as.integer(N),
               as.integer(iter),
               PACKAGE = "dpOMP")

  names(out) <- c("beta", "pi", "beta_g", "sigma2")
  out$beta <- array(out$beta, dim=c(V, K, iter))
  out$pi <- array(out$pi, dim=c(K, iter))
  out$beta_g <- array(out$beta_g, dim=c(V, G, iter))
  return(out)
}