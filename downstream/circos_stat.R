library(tidyverse)


CG <- read_delim('../data/circos_data/WNYCY.CG.bed',delim = '\t',col_names = F)
CHG <- read_delim('../data/circos_data/WNYCY.CHG.bed',delim = '\t',col_names = F)
CHH <- read_delim('../data/circos_data/WNYCY.CHH.bed',delim = '\t',col_names = F)
TE <- read_delim('../data/circos_data/WNYCY.EDTA.TE.bed',delim = '\t',col_names = F)
gene <- read_delim('../data/circos_data/WNYCY.gene.high_confidence.bed',delim = '\t',col_names = F)
CENH3 <- read_delim('../data/circos_data/WNYCY.CENH3.bed',delim = '\t',col_names = F)


ggplot(CG) + geom_density(aes(x = X4))  + xlim(c(0.8,1))
ggplot(CHG) + geom_density(aes(x = X4)) + xlim(c(0.4,0.8))
ggplot(CHH) + geom_density(aes(x = X4)) + xlim(c(0,0.05))
ggplot(TE) + geom_density(aes(x = X4)) + xlim(c(0,200))
ggplot(gene) + geom_density(aes(x = X4)) + xlim(c(0,30))
ggplot(CENH3) + geom_density(aes(x = X4)) + xlim(c(0,0.4))
