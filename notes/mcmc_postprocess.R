source("plot_mcmc.R")

mcmc <- readRDS("july05_1d/samples.rds")
data <- readRDS("july05_1d/data.rds")
make_ppplot_for_P(mcmc = mcmc, true = data$beta, dims=data$V, names="beta", plotfilename = "ppplot_testrun")