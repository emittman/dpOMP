source("plot_mcmc.R")

mcmc <- readRDS("july05_1d/samples_with_inits.rds")
data <- readRDS("july05_1d/data.rds")
make_ppplot_for_P(mcmc = mcmc, true = data$beta, dims=1, names="beta", plotfilename = "ppplot_testrun")