# Functions for summaries and plots

#is this actually posterior predictive?
make_ppplot_for_P <- function(mcmc, true.beta = NULL, true.pi = NULL, dims, names, bins=30, plotfilename){
  require(ggplot2)
  #mcmc$beta_g is V by K by n_iter array
  if(dims != length(names))
    error("Length of dims much match length of names")
  
  if(dims < 1)
    error("dims must be at least 1")
    
  if(dims == 1){
    plot.df <- data.frame(x1 = as.numeric(mcmc$beta_g))
    p <- ggplot(plot.df, aes(x=x1)) + geom_histogram(aes(y=..density..), bins=bins) +
      xlab(names)
    # if(!is.null(true.beta)){
      ref.df <- data.frame(x1 = as.numeric(true.beta), pi = as.numeric(true.pi))
      p <- p + geom_segment(data=ref.df, aes(x = x1, xend = x1, y=0, yend=pi), color = "red")
    # }
    pdf(paste0(c(plotfilename, ".pdf"), collapse=""))
    print(p)
    dev.off()
  }
  
  if(dims > 1){
    dims <- dims[1:2]
    names <- names[1:2]
    require(hexbin)
    plot.df <- data.frame(x1 = as.numeric(mcmc$beta_g[1, , ]), x2 = as.numeric(mcmc$beta_g[2, , ]))
    p <- ggplot(plot.df, aes(x=x1, y=x2)) + geom_hex(bins=bins) +
      xlab(names[1]) + ylab(names[2]) +
      scale_fill_continuous(trans="log")
    # if(!is.null(true)){
      ref.df <- data.frame(x1 = as.numeric(true.beta[1, ]), x2 = as.numeric(true.beta[2, ]), pi = as.numeric(true.pi))
      p <- p + geom_point(data = ref.df, aes(x = x1, y = x2, size=pi), color = "red")
    # }
    pdf(paste0(c(plotfilename, ".pdf"), collapse=""))
    print(p)
    dev.off()
  }
}

