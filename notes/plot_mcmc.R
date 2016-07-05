# Functions for summaries and plots

#is this actually posterior predictive?
make_ppplot_for_P <- function(mcmc, dims, names, bins=30, plotfilename){
  require(ggplot2)
  #mcmc$beta_g is V by K by n_iter array
  len = length(dims)
  if(len != length(names))
    error("Length of dims much match length of names")
  
  if(len < 1)
    error("dims must be at least 1")
    
  if(len == 1){
    plot.df <- data.frame(x1 = as.numeric(mcmc$beta_g))
    p <- ggplot(plot.df, aes(x=x1)) + geom_histogram(bins=bins) +
      xlab(name)
    print(p)
    ggsave(paste0(c(plotfilename, ".pdf"), collapse=""))
  }
  
  if(len > 1){
    dims <- dims[1:2]
    names <- names[1:2]
    require(hexbin)
    plot.df <- data.frame(x1 = as.numeric(mcmc$beta_g[dim[1], , ]), x2 = as.numeric(mcmc$beta_g[dim[2], , ]))
    p <- ggplot(plot.df, aes(x=x1, y=x2)) + geom_hex(bins=bins) +
      xlab(name[1]) + ylab(name[2]) +
      scale_fill_continuous(trans="log")
    print(p)
    ggsave(paste0(c(plotfilename, ".pdf"), collapse=""))
  }
}

