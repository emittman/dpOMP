source("plot_mcmc.R")

mcmc <- readRDS("samples_with_inits.rds")
make_ppplot_for_P(mcmc = mcmc, dims=1, names="beta", plotfilename = "ppplot_testrun")