library(tidyverse)
library(coda)
library(GGally)

out <- readRDS("Outputs/VegSeasonYear_TE_100K.rds")
samples <- out$samples

## pairs plot
samples <- as.matrix(samples)
samples <- as.data.frame(samples)
names(samples)
samples <- samples %>%
  select("available_prob[1]", "available_prob[2]", "beta[1]", "beta[2]", "beta[3]", "beta[4]",
         "beta[5]","beta[6]","beta[7]","beta[8]", "beta[9]", "beta[10]", "beta[11]", "sigma[1]", "sigma[2]")

colnames(samples) <- c("Pm", "a1", "a2", "b1", "b2", "betaa1", "betaa2", "betac1", "c1")
samples <- samples[sample.int(nrow(samples), ceiling(nrow(samples) * 0.1)), ]
samples %>%
  as.data.frame() %>%
  ggpairs()
