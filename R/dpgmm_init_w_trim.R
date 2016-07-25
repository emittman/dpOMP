#' Run MCMC for DPGMMM model, utilizing pilot chain only saving samples for certain genes
#' 
#' @param data A column-major $\eqn{G \cdot VN}$ gene expressionvector/matrix
#' @param design A column-major $\eqn{V \cdot V}$ experimental design vector/matrix
#' @param lambda2 prior variance for base measure
#' @param alpha DP concentration parameter
#' @param G number of genes
#' @param K maximum number of clusters
#' @param V number of varieties/groups
#' @param N samples per variety/group
#' @param iter number of iterations for main chain
#' @param init_iter number of iterations for pilot chain
#' @param save_index zero-indexed vector indicating for which genes the save parameters will be saved
#' 
#' @export
dpgmm_init_w_trim <- function(data, design, lambda2, alpha, a, b,  G, V, K, N, iter, init_iter, save_index){
  # run an initial chain
  init_out <- dpgmm(data=data, design=design, lambda2=lambda2, alpha=alpha,
                    a=a, b=b, G=G, V=V, K=K, N=N, iter=init_iter)
  
  yTy <- apply(data, 1, crossprod)
  #yTy <- sum(yTy)
  xTy <- t(design) %*% t(data)
  xTx <- crossprod(design)
  
  # move most relevant atoms to front of vector
  init <- list()
  indices <- which(init_out$beta[1,,init_iter] %in% init_out$beta_g[1,,init_iter])
  anti_indices <- (1:K)[!(1:K %in% indices)]
  indices <- indices[order(-init_out$pi[indices,init_iter])]
  init$beta <- init_out$beta[,c(indices, anti_indices),init_iter]
  init$sigma2 <- init_out$sigma2[c(indices, anti_indices),init_iter]
  init$pi <- init_out$pi[c(indices, anti_indices),init_iter]
  
  num_save <- length(save_index)
  
  out <- .Call("dpgmm_init_w_trimR",
               as.numeric(yTy),
               as.numeric(xTy),
               as.numeric(xTx),
               as.numeric(init$beta),
               as.numeric(init$pi),
               as.numeric(init$sigma2),
               as.numeric(lambda2),
               as.numeric(alpha),
               as.numeric(a),
               as.numeric(b),
               as.integer(G),
               as.integer(V),
               as.integer(K),
               as.integer(N),
               as.integer(iter),
               as.integer(save_index),
               as.integer(num_save),
               PACKAGE = "dpOMP")
  
  names(out) <- c("beta_g", "sigma2_g", "max_index")
  #   out$beta <- array(out$beta, dim=c(V, K, iter))
  #   out$pi <- array(out$pi, dim=c(K, iter))
  #   out$sigma2 <- array(out$sigma2, dim=c(K, iter))
  out$beta_g <- array(out$beta_g, dim=c(V, num_save, iter))
  out$sigma2_g <- array(out$sigma2_g, dim=c(num_save, iter))
  return(out)
}