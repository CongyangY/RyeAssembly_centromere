library(tidyverse)
library(reshape2)
library(RColorBrewer)
library(RIdeogram)

# 0 preparation -----------------------------------------------------------

library(showtext)
showtext_auto()
# 如没有helvetica字体，可以考虑使用Califri，二者很相似


# windowsFonts(helvetica=windowsFont("Helvetica CE 55 Roman"),
#              Times=windowsFont("Times New Roman"),
#              Arial=windowsFont("Arial"))
# 将windows字体映射到ggplot中
# 

theme_Publication <- function(base_size=14, base_family="helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(size = 12), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(0.2, "cm"),
            legend.margin = margin(unit(0, "cm")),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}

scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}



# 也可以使用下面两个函数
scale_fill_Publication <- function(color_name = 'Set2'){
  library(scales)
  library(dplyr)
  library(RColorBrewer)
  test_color <- as.data.frame(brewer.pal.info)
  test_color$name <- rownames(test_color)
  special_color <- filter(test_color,name == color_name)
  discrete_scale("fill","Publication",manual_pal(values = brewer.pal(special_color[,1],special_color[,4])))
  
}

scale_colour_Publication <- function(color_name = 'Set2'){
  library(scales)
  library(RColorBrewer)
  library(RColorBrewer)
  test_color <- as.data.frame(brewer.pal.info)
  test_color$name <- rownames(test_color)
  special_color <- filter(test_color,name == color_name)
  discrete_scale("colour","Publication",manual_pal(values = brewer.pal(special_color[,1],special_color[,4])))
  
}






# centromere location -----------------------------------------------------
WN_karyotype <- 
  data.frame(
    Chr = c(paste('chr',1:7,'R',sep = '')),
    Start = c(rep(0,7)),
    End = c(991023136,1118204079,1156965858,1093317930,
            1257695468,1165049686,999733186),
    CE_start = c(429700000,474000000,537000000,387000000,460000000,
                 366000000,471000000),
    CE_end = c(440600000,483000000,549500000,398600000,471500000,
               375800000,488000000)
  )

ideogram(karyotype = WN_karyotype)
file.rename('chromosome.svg','../figure/CENH3/WN.chromosome.svg')
convertSVG("../figure/CENH3/WN.chromosome.svg", device = "png")

# CENH3-enriched/depleted region -------------------------------------------------------------


CENH3 <- read_delim('../data/CENH3/wn_CENH3.cent.bg',delim = '\t',col_names = F)

Input <- read_delim('../data/CENH3/wn_input.bg',delim = '\t',col_names = F)
mean(Input$X4) # 0.1342459
sd(Input$X4) # 0.08946951

CENH3$p <- apply(CENH3[,4], 1, function(x){pnorm(x,mean = 0.1342459,sd = 0.08946951,lower.tail = F)})



write_delim(filter(CENH3,p<=1e-5)[,1:3],'../data/CENH3/CENH3_enriched.bed',delim = '\t',col_names = F)
write_delim(filter(CENH3,p>1e-5)[,1:3],'../data/CENH3/CENH3_depleted.bed',delim = '\t',col_names = F)
write_delim(CENH3,'../data/CENH3/CENH3_enrich_pvalue.bed',delim = '\t',col_names = F)



CENH3_enrich <- CENH3 %>% filter(p<=1e-5) %>% group_by(X1) %>% 
  mutate(diff = X3-X2) %>%
  summarise(Enrich = sum(diff))
CENH3_deplete <- CENH3 %>% filter(p>1e-5) %>% group_by(X1) %>% 
  mutate(diff = X3-X2) %>%
  summarise(Deplete = sum(diff))
CENH3_subdomain <- full_join(CENH3_enrich,CENH3_deplete,by = join_by(X1 == X1)) %>% 
  melt(id = 'X1')
ggplot(CENH3_subdomain) + geom_bar(aes(x = X1, y = value/1e6,fill = variable,color = variable),
                                   stat = 'identity',position = 'stack',width = 0.6) + 
  scale_colour_Publication() + scale_fill_Publication() + 
  theme_Publication() + labs(
    x = 'Chromosome',
    y = 'Length (Mb)',
    fill = 'Subdomain',
    color = 'Subdomain'
  ) + 
  theme(
    legend.position = 'top'
  )
ggsave('../figure/CENH3/CENH3.subdomain.bar.pdf',width = 5,height = 4)




# centromere sequence type ------------------------------------------------

seq.type <- read_delim('../data/CENH3/WNYCY.cent.seq_type.list',delim = '\t',col_names = T) %>% 
  melt(id = 'Chr')
ggplot(seq.type) + geom_bar(aes(x = Chr, y = value/1e6, fill = variable,color = variable),
                            position = 'stack',stat = 'identity',width = 0.6) + 
  theme_Publication() + scale_fill_Publication('Set2') + scale_colour_Publication('Set2') + 
  labs(x = 'Chromosome', y = 'Length (Mb)', fill = 'Type',color = 'Type')
ggsave('../figure/CENH3/CENH3.LTR.length.bar.pdf',width = 4,height = 4)




# CENH3 signal  -----------------------------------------------------------

cenh3 <- read_delim('../data/CENH3/wn_CENH3.100k.cent.bg',delim = '\t',col_names = F) %>% 
  mutate(Type = 'CENH3')
input <- read_delim('../data/CENH3/wn_input.100k.cent.bg',delim = '\t',col_names = F) %>% 
  mutate(Type = 'Input')

signal <- rbind(cenh3,input)

ggplot(filter(signal,X1 == 'chr1R')) + 
  geom_line(aes(x = (X3+X2)/2,y=X4,color =  Type)) + 
  theme_Publication() + 
  labs(x = 'Position', y = 'Mean signal level') + 
  theme(
    legend.position = 'none',
    legend.direction = 'vertical',
    axis.text = element_blank(),
    axis.ticks = element_blank(),
  ) + labs(x = '', y = '')
ggsave('../figure/CENH3/chr1.cenh3_signal.meta.pdf',width = 30,height = 3)

ggplot(filter(signal,X1 == 'chr2R')) + 
  geom_line(aes(x = (X3+X2)/2,y=X4,color =  Type)) + 
  theme_Publication() + 
  labs(x = 'Position', y = 'Mean signal level') + 
  theme(
    legend.position = 'none',
    legend.direction = 'vertical',
    axis.text = element_blank(),
    axis.ticks = element_blank(),
  ) + labs(x = '', y = '')
ggsave('../figure/CENH3/chr2.cenh3_signal.meta.pdf',width = 30,height = 3)

ggplot(filter(signal,X1 == 'chr3R')) + 
  geom_line(aes(x = (X3+X2)/2,y=X4,color =  Type)) + 
  theme_Publication() + 
  labs(x = 'Position', y = 'Mean signal level') + 
  theme(
    legend.position = 'none',
    legend.direction = 'vertical',
    axis.text = element_blank(),
    axis.ticks = element_blank(),
  ) + labs(x = '', y = '')
ggsave('../figure/CENH3/chr3.cenh3_signal.meta.pdf',width = 30,height = 3)

ggplot(filter(signal,X1 == 'chr4R')) + 
  geom_line(aes(x = (X3+X2)/2,y=X4,color =  Type)) + 
  theme_Publication() + 
  labs(x = 'Position', y = 'Mean signal level') + 
  theme(
    legend.position = 'none',
    legend.direction = 'vertical',
    axis.text = element_blank(),
    axis.ticks = element_blank(),
  ) + labs(x = '', y = '')
ggsave('../figure/CENH3/chr4.cenh3_signal.meta.pdf',width = 30,height = 3)

ggplot(filter(signal,X1 == 'chr5R')) + 
  geom_line(aes(x = (X3+X2)/2,y=X4,color =  Type)) + 
  theme_Publication() + 
  labs(x = 'Position', y = 'Mean signal level') + 
  theme(
    legend.position = 'none',
    legend.direction = 'vertical',
    axis.text = element_blank(),
    axis.ticks = element_blank(),
  ) + labs(x = '', y = '')
ggsave('../figure/CENH3/chr5.cenh3_signal.meta.pdf',width = 30,height = 3)

ggplot(filter(signal,X1 == 'chr6R')) + 
  geom_line(aes(x = (X3+X2)/2,y=X4,color =  Type)) + 
  theme_Publication() + 
  labs(x = 'Position', y = 'Mean signal level') + 
  theme(
    legend.position = 'none',
    legend.direction = 'vertical',
    axis.text = element_blank(),
    axis.ticks = element_blank(),
  ) + labs(x = '', y = '')
ggsave('../figure/CENH3/chr6.cenh3_signal.meta.pdf',width = 30,height = 3)

ggplot(filter(signal,X1 == 'chr7R')) + 
  geom_line(aes(x = (X3+X2)/2,y=X4,color =  Type) ) + 
  theme_Publication() + 
  labs(x = 'Position', y = 'Mean signal level') + 
  theme(
    legend.position = 'none',
    legend.direction = 'vertical',
    axis.text = element_blank(),
    axis.ticks = element_blank(),
  ) + labs(x = '', y = '')
ggsave('../figure/CENH3/chr7.cenh3_signal.meta.pdf',width = 30,height = 3)





