library(Paschold2012)
library(plyr)
library(dplyr)
library(dpOMP)
library(ggplot2)

#random subset of genes
G <- 10000
genes <- sample(unique(Paschold2012$GeneID), size=G, replace=FALSE)

data <- Paschold2012 %>% 
  subset(genotype %in% c("B73", "Mo17") & GeneID %in% genes) %>%
  arrange(GeneID, genotype, replicate)

head(data, 20)
summary(data)



V <- 2
N <- 4
X <- matrix(rep(c(1, 1, -1, 1), each=4), V*N, V)
y <- matrix(log(data$total+1), byrow=T, G, V*N)

out <- dpgmm_init(data=y, design=X, G=G, V=V, K=100, N=N, iter=10000, init_iter=3000)
saveRDS(out, file="diffexpr_628.rds")

out <- readRDS("diffexpr_628.rds")
#number of occupied clusters
occupied <- sapply(1:10000, function(i) sum(out$beta[1,,i] %in% out$beta_g[1,,i]))
hist(occupied, breaks = 1:100+.5)

#maximum index of significant pi
max_index <- sapply(1:10000, function(i) max(which(out$pi[,i]>.005)))
hist(max_index, breaks = 1:100+.5)

#plot beta_g's
Bhat <- data.frame(V1 = as.numeric(out$beta_g[1,,10000]), V2 = as.numeric(out$beta_g[2,,9901:10000]))

ggplot(Bhat, aes(x=V1, y=V2)) + geom_hex(bins=30)
