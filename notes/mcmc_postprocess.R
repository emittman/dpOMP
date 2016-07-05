source("plot_mcmc.R")

mcmc <- readRDS("samples_with_inits.rds")
make_ppplot_for_P(mcmc = mcmc, dim1 = 1, dim2 = 2, name1 = "X1", name2 = "X2", plotfilename = "ppplot_testrun")