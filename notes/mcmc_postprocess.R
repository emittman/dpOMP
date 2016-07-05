source("plot_mcmc.R")

mcmc <- readRDS("july05_1d/samples.rds")
data <- readRDS("july05_1d/data.rds")
make_ppplot_for_P(mcmc = mcmc, true.beta = data$beta, true.pi = data$pi, dims=data$V, names="beta", plotfilename = "july05_1d/ppplot_testrun")

mcmc <- readRDS("july05_2d/samples.rds")
data <- readRDS("july05_2d/data.rds")
make_ppplot_for_P(mcmc = mcmc, true.beta = data$beta, true.pi = data$pi, dims=data$V, names=c("beta1","beta2"), plotfilename = "july05_2d/ppplot_testrun")
