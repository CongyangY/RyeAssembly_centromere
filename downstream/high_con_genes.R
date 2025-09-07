library(tidyverse)

featureC <- read_delim('../data/WNhan.EVM.featureCount',delim = '\t',skip = 2,col_names = F)

featureC[,7:68] <- 10^6*featureC[,7:68] /colSums(featureC[,7:68])
featureC <- featureC %>% mutate(
  across(7:68,~./featureC[[6]]*10^3) # cross(7:68,~./featureC[,6]*10^3),error
) %>% filter(rowSums(select(., 8:63) > 0.1) > 0) 


tmp <- featureC%>% filter(rowSums(select(., 8:63) > 0.1) > 5) 

featureC[,7:68] <- as.numeric(featureC[,7:68])
ggplot(featureC) + 
  geom_density(aes(x = log10(X8)))
