library(tidyverse)
library(reshape2)
library(RColorBrewer)


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
# length distribution of tandem repeat  -----------------------------------


Set1 <- brewer.pal(8,'Set1')

tr.571 <- read_delim('../data/repeat_consensus/repeatExplorer.bed',delim = '\t',col_names = F) %>% 
  filter(X4=='repeat_571',X3-X2>200) %>% 
  mutate(
    len = X3-X2,
    color = ifelse(len == 571,'max','other'))

ggplot(tr.571) + geom_bar(aes(x = len,color = color,fill = color),width = 0.5) + xlim(x = c(500,650)) + 
  theme_Publication() + scale_color_manual(values = c('max' = Set1[1],'other' = 'gray')) +
  scale_fill_manual(values = c('max' = Set1[1],'other' = 'gray')) + 
  labs(x = 'Length',y = 'COunt') + theme(legend.position = 'None') + 
  geom_bar(data = filter(tr.571,color !='other'),width = 0.3 ,aes(x = len,color = color,fill = color))
ggsave('../figure/repeat_consensus/tr571.len.dis.pdf',width = 4.6,height = 3.3)


tr.442 <- read_delim('../data/repeat_consensus/pSc.bed',delim = '\t',col_names = F) %>% 
  filter(X4=='pSc119.2',X3-X2>200) %>% 
  mutate(
    len = X3-X2,
    color = ifelse(len %in% c(413,442,439),'max','other'))

ggplot(tr.442) + geom_bar(aes(x = len,color = color,fill = color),width = 0.5) + xlim(x = c(350,500)) + 
  theme_Publication() + scale_color_manual(values = c('max' = Set1[1],'other' = 'gray')) +
  scale_fill_manual(values = c('max' = Set1[1],'other' = 'gray')) + 
  labs(x = 'Length',y = 'COunt') + theme(legend.position = 'None') + 
  geom_bar(data = filter(tr.442,color !='other'),width = 0.3 ,aes(x = len,color = color,fill = color))
ggsave('../figure/repeat_consensus/tr442.len.dis.pdf',width = 4.6,height = 3.3)




tr.520 <- read_delim('../data/repeat_consensus/pSc.bed',delim = '\t',col_names = F) %>% 
  filter(X4=='pSc200',X3-X2>200) %>% 
  mutate(
    len = X3-X2,
    color = ifelse(len %in% c(511,520),'max','other'))

ggplot(tr.520) + geom_bar(aes(x = len,color = color,fill = color),width = 0.3) + xlim(x = c(400,600)) + 
  theme_Publication() + scale_color_manual(values = c('max' = Set1[1],'other' = 'gray')) +
  scale_fill_manual(values = c('max' = Set1[1],'other' = 'gray')) + 
  labs(x = 'Length',y = 'COunt') + theme(legend.position = 'None') + 
  geom_bar(data = filter(tr.520,color !='other'),width = 0.3 ,aes(x = len,color = color,fill = color))
ggsave('../figure/repeat_consensus/tr.520.len.dis.pdf',width = 4.6,height = 3.3)



tr.571 <- read_delim('./tmp1.cons.bed',delim = '\t',col_names = F) %>% 
  mutate(
    len = X3-X2,
    color = ifelse(len == 572,'max','other'))

ggplot(tr.571) + geom_bar(aes(x = len,color = color,fill = color),width = 0.5) + xlim(x = c(500,650)) + 
  theme_Publication() + scale_color_manual(values = c('max' = Set1[1],'other' = 'gray')) +
  scale_fill_manual(values = c('max' = Set1[1],'other' = 'gray')) + 
  labs(x = 'Length',y = 'COunt') + theme(legend.position = 'None') + 
  geom_bar(data = filter(tr.571,color !='other'),width = 0.3 ,aes(x = len,color = color,fill = color))
