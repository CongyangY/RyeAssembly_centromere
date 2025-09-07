library(tidyverse)
library(reshape2)

# 0 preparation -----------------------------------------------------------

library(showtext)

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




aver <- function(data){mean(data,na.rm = T,trim = 0.0025)}
Dark2 <- brewer.pal(8,'Dark2')
Paired <- brewer.pal(8,'Paired')
Set2 <- brewer.pal(8,'Set2')
Set1 <- brewer.pal(8,'Set1')
# path: 读取的文件路径
# names： 每个line所读应的name，要与bw文件顺序一一对应
# n_max: 默认读取全部行。当文件过大时，可以自定义读取的行数
plot_data <- function(path,names,n_max=Inf){
  histone <- read_delim(path, delim = '\t', skip =1, n_max = n_max,col_names = F)
  
  nums = length(names) # bw文件的个数
  data <- as.data.frame(matrix(rep(NA,dim(histone)[[2]]-6),ncol = nums))
  for (i in 1:nums){
    # total bin  
    sum_bins = dim(histone)[[2]] - 6
    each_bins = sum_bins/nums
    
    start = 7 + (i-1)*each_bins
    end = start + each_bins -1
    data[,i] <- apply(histone[,start:end],2,aver)
  }
  
  
  pos <- 1:each_bins
  data <- cbind(data,pos)
  colnames(data) <- c(names,'pos')
  
  result <- melt(data,id=(c('pos')))
  
  return(result)
}



# gene length -------------------------------------------------------------


arm.gene.len <- read_delim('../data/gene/WNhan.EVM.high_confidence.noN.ArmGene.length',
                           delim = '\t',col_names = F) %>% mutate(Type = 'Arm')
cent.gene.len <- read_delim('../data/gene/WNhan.EVM.high_confidence.noN.centGene.length',
                           delim = '\t',col_names = F) %>% mutate(Type = 'Cent')

gene.len <- rbind(arm.gene.len,cent.gene.len)
ggplot(gene.len) + geom_density(aes(x = X1,color = Type)) + 
  theme_Publication() + scale_colour_Publication('Dark2') + 
  labs(y = 'Density', x = 'Length(bp)') + theme(
    legend.position = c(0.7,0.7),
    legend.direction = 'vertical'
  )
ggsave('../figure/gene/GeneLen.density.pdf',width = 4,height = 2.8)

ggplot(gene.len,aes(y = log10(X1),x = Type,fill = Type)) + geom_boxplot() + 
  geom_signif(comparisons = list(c("Arm", "Cent")), y_position = 4.8,map_signif_level=F) + 
  theme_Publication() + labs(x = 'Type',y = 'log10(Length (bp))') + 
  theme(
    legend.position = 'None'
  )
ggsave('../figure/gene/GeneLen.boxplot.pdf',width = 2.5,height = 2.8)



# centromere region  ------------------------------------------------------

cenh.region <- data.frame(
  Type = c('Enriched','Depleted'),
  Nums = c(10,15)
)

ggplot(cenh.region) + geom_bar(aes(x = Type, y = Nums,fill = Type),stat = 'identity')


cenh3.enrich <- c(1,2,3,9,10,15,21,22,23,24)
cenh3.deplete <- c(4,5,6,7,8,11,12,13,14,16,17,18,19,20,25)
inactive.gene <- c(4,5,6,12,17,18,19)
result = list(
  'Enriched-region' = cenh3.enrich,
  'Depleted-region' = cenh3.deplete,
  'Activated gene' = inactive.gene
)

library('ggvenn')


ggvenn(result)
ggsave('../figure/gene/centRegionGene.inactive.venn.pdf',
       width = 4.5,height = 3.5)
