# Functions for summaries and plots

#is this actually posterior predictive?
make_ppplot_for_P <- function(mcmc, dim1, dim2, name1, name2, bins, plotfilename){
  #mcmc$beta_g is V by K by n_iter array
  require(ggplot2)
  require(hexbin)
  plot.df <- data.frame(x1 = mcmc$beta_g[dim1, , ], x2 = mcmc$beta_g[dim2, , ])
  p <- ggplot(plot.df, aes(x=x1, y=x2)) + geom_hex(bins=bins) +
    xlab(name1) + ylab(name2) +
    scale_fill_continuous(trans="log")
  print(p)
  ggsave(paste0(c(plotfilename, ".pdf"), collapse=""))
}

