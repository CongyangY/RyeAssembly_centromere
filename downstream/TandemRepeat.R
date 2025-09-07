library(tidyverse)
library(RColorBrewer)

# 0. preparation ----------------------------------------------------------

library(showtext)

# 如没有helvetica字体，可以考虑使用Califri，二者很相似


# windowsFonts(helvetica=windowsFont("Helvetica CE 55 Roman"),
#              Times=windowsFont("Times New Roman"),
#              Arial=windowsFont("Arial"))
# 将windows字体映射到ggplot中
# 

theme_Publication <- function(base_size=14) {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size)
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




# 1. length distribution----------------------------------------------------------------------
Set1 <- brewer.pal(8,'Set1')


TR_442 <- read_delim('../data/TandemRepeat/442.bed',col_names = F) %>% 
  mutate(len = X3-X2,color = ifelse(len == 442,'max','other'))
ggplot(TR_442) + geom_bar(aes(x = len,color = color,fill = color),width = 0.5)  + 
  theme_Publication() + scale_color_manual(values = c('max' = Set1[1],'other' = 'gray')) +
  scale_fill_manual(values = c('max' = Set1[1],'other' = 'gray')) + 
  labs(x = 'Length',y = 'Count') + theme(legend.position = 'None') + 
  geom_bar(data = filter(TR_442,color !='other'),width = 0.3 ,aes(x = len,color = color,fill = color))
ggsave('../figure/TandemRepeat/TR_442.len.dis.pdf',width = 3.6,height = 2)


TR_520 <- read_delim('../data/TandemRepeat/520.bed',col_names = F) %>% 
  mutate(len = X3-X2,color = ifelse(len == 520,'max','other'))
ggplot(TR_520) + geom_bar(aes(x = len,color = color,fill = color),width = 0.5) + xlim(c(450,550))  + 
  theme_Publication() + scale_color_manual(values = c('max' = Set1[1],'other' = 'gray')) +
  scale_fill_manual(values = c('max' = Set1[1],'other' = 'gray')) + 
  labs(x = 'Length',y = 'COunt') + theme(legend.position = 'None') + 
  geom_bar(data = filter(TR_520,color !='other'),width = 0.3 ,aes(x = len,color = color,fill = color))
ggsave('../figure/TandemRepeat/TR_520.len.dis.pdf',width = 3.6,height = 2)


TR_571 <- read_delim('../data/TandemRepeat/571.bed',col_names = F) %>% 
  mutate(len = X3-X2,color = ifelse(len == 571,'max','other'))
ggplot(TR_571) + geom_bar(aes(x = len,color = color,fill = color),width = 0.5) + xlim(c(500,620))  + 
  theme_Publication() + scale_color_manual(values = c('max' = Set1[1],'other' = 'gray')) +
  scale_fill_manual(values = c('max' = Set1[1],'other' = 'gray')) + 
  labs(x = 'Length',y = 'COunt') + theme(legend.position = 'None') + 
  geom_bar(data = filter(TR_571,color !='other'),width = 0.3 ,aes(x = len,color = color,fill = color))
ggsave('../figure/TandemRepeat/TR_571.len.dis.pdf',width = 3.6,height = 3.3)


TR_380 <- read_delim('../data/TandemRepeat/380.all.bed',col_names = F) %>% 
  mutate(len = X3-X2,color = ifelse(len == 380,'max','other'))
ggplot(TR_380) + geom_bar(aes(x = len,color = color,fill = color),width = 0.5) + xlim(0,420) + 
  theme_Publication() + scale_color_manual(values = c('max' = Set1[1],'other' = 'gray')) +
  scale_fill_manual(values = c('max' = Set1[1],'other' = 'gray')) + 
  labs(x = 'Length',y = 'Count') + theme(legend.position = 'None') + 
  geom_bar(data = filter(TR_380,color !='other'),width = 0.3 ,aes(x = len,color = color,fill = color))
ggsave('../figure/TandemRepeat/TR_380.len.dis.pdf',width = 3.6,height = 3.3)


# 2. variation ------------------------------------------------------------

pos_cov.442 <- read_delim('../data/TandemRepeat/442_cons.var',
                          delim = '\t',col_names = F)
pos_cov.520 <- read_delim('../data/TandemRepeat/520_cons.var',
                          delim = '\t',col_names = F)
pos_cov.571 <- read_delim('../data/TandemRepeat/571_cons.var',
                          delim = '\t',col_names = F)
pos_cov.380 <- read_delim('../data/TandemRepeat/380.1k.var',
                          delim = '\t',col_names = F)
ggplot(pos_cov.442) + geom_bar(aes(x = X1,y=X2),stat = 'identity') + ylim(c(0,0.8)) + 
  theme_Publication() + labs(x = 'Position (bp)', y = 'Variant frequency')
ggsave('../figure/TandemRepeat/TR_442.pos_cov.barplot.pdf',width = 3.6, height = 2)

ggplot(pos_cov.520) + geom_bar(aes(x = X1,y=X2),stat = 'identity') + ylim(c(0,0.8)) + 
  theme_Publication() + labs(x = 'Position (bp)', y = 'Variant frequency')
ggsave('../figure/TandemRepeat/TR_520.pos_cov.barplot.pdf',width = 3.6, height = 2)

ggplot(pos_cov.571) + geom_bar(aes(x = X1,y=X2),stat = 'identity') + ylim(c(0,0.8)) + 
  theme_Publication() + labs(x = 'Position (bp)', y = 'Variant frequency')
ggsave('../figure/TandemRepeat/TR_571.pos_cov.barplot.pdf',width = 3.6, height = 2)

ggplot(pos_cov.380) + geom_bar(aes(x = X1,y=X2),stat = 'identity') + ylim(c(0,0.8)) + 
  theme_Publication() + labs(x = 'Position (bp)', y = 'Variant frequency')
ggsave('../figure/TandemRepeat/TR_380.pos_cov.barplot.pdf',width = 3.6, height = 2)


# degrade -----------------------------------------------------------------

test <- read_delim('/Users/ycy/Downloads/all_442.2.matrix',delim = '\t',col_names = F)
library(pheatmap)
test <- as.matrix(test)
pheatmap(test, cluster_row = F,cluster_col = FALSE,show_colnames = F,
         show_rownames = F,width=3,height=5,scale = 'row',angle_col = "45",
         breaks = seq(0.5, 1, length.out = 100))
test2 <- read_delim('/Users/ycy/Downloads/all_520.2.matrix',delim = '\t',col_names = F)
test2 <- as.matrix(test2)
pheatmap(test2, cluster_row = F,cluster_col = FALSE,show_colnames = F,
         show_rownames = F,width=3,height=5,scale = 'row',angle_col = "45",
         breaks = seq(0.5, 1, length.out = 100))


test3 <- read_delim('/Users/ycy/Downloads/out_matrix.csv',delim = ',',col_names = T)
test3 <- test3[,-1]
test3 <- as.matrix(test3)
pheatmap(test3, cluster_row = F,cluster_col = FALSE,show_colnames = F,
         show_rownames = F,width=3,height=5,scale = 'row',angle_col = "45",
         breaks = seq(0.5, 1, length.out = 100))

test4 <- read_delim('/Users/ycy/Downloads/571.filter_plot.matrix',delim = '\t',col_names = T)
test4 <- test4[,-1]
test4 <- as.matrix(test4)
pheatmap(test4, cluster_row = F,cluster_col = FALSE,show_colnames = F,
         show_rownames = F,width=3,height=5,scale = 'row',angle_col = "45",
         breaks = seq(0.8, 1, length.out = 100))



# LTR density -------------------------------------------------------------

TR <- read_delim('../data/TandemRepeat/WN.2M.TR_LTR.bed',delim = '\t',col_names = F) %>% 
  mutate(Type = 'Subtelemere')
Arm <- read_delim('../data/TandemRepeat/WN.2M.ARM_LTR.bed',delim = '\t',col_names = F) %>% 
  mutate(Type = 'Arm')
Peri <- read_delim('../data/TandemRepeat/WN.2M.Peri_LTR.bed',delim = '\t',col_names = F) %>% 
  mutate(Type = 'Peri')
Cent <- read_delim('../data/TandemRepeat/WN.2M.Cent_LTR.bed',delim = '\t',col_names = F) %>% 
  mutate(Type = 'Cent')

LTR <- rbind(TR,Arm,Peri,Cent)
LTR$Type <- factor(LTR$Type,levels = c('Subtelemere','Arm','Peri','Cent'))


library('ggsignif')

ggplot(LTR,aes(x = Type, y = X4)) + geom_boxplot(aes(fill = Type),outlier.size = 0.1) + 
  theme_Publication() + scale_fill_Publication('Set1') + 
  geom_signif(comparisons=list(c("Subtelemere", "Arm")),map_signif_level = T,y_position = 100) + 
  geom_signif(comparisons=list(c("Peri", "Arm")),map_signif_level = T,y_position = 110) + 
  geom_signif(comparisons=list(c("Peri", "Cent")),map_signif_level = T,y_position = 130) + 
  theme(legend.position = 'None') + labs(
    x = 'Region', y = 'Relative Frequency of full-length LTR'
  )
ggsave('../figure/Region_LTR.density.boxplot.pdf',width = 5,height = 6)
  

LTR <- read_delim('../data/TandemRepeat/WN.2M.TR_LTR.bed',delim = '\t',col_names = F) %>% 
  mutate(Type = 'LTR')
H3K27 <- read_delim('../data/TandemRepeat/WN.2M.TR_H3K27me3.bed',delim = '\t',col_names = F) %>% 
  mutate(Type = 'H3K27me3',Groyp = ifelse(X4<=0.1,'Low',ifelse(X4<=0.25,))) 
H3K9 <- read_delim('../data/TandemRepeat/WN.2M.TR_H3K9me3.bed',delim = '\t',col_names = F) %>% 
  mutate(Type = 'H3K9me3')

LTR_K9 <- data.frame(
  LTR = LTR$X4,
  K9 =H3K9$X4
)
summary(lm(log10(H3K9$X4)~log10(LTR$X4+0.1)))
ggplot(LTR_K9) + geom_abline(slope = 0.67092,intercept = -1.92543,color = Set1[1],size = 1) +
  geom_jitter(aes(x = log10(LTR),y = log10(K9))) + 
  theme_Publication() + 
  labs(x = 'LTR density (log10)',y = 'H3K9me3 intensity (log10)') +
  annotate("text", x = 0.05, y = max(log10(LTR_K9$K9))+0.5, 
           label = paste("y = 0.67x - 1.93\nR² = 0.21, P < 0.001"),
           hjust = 0, vjust = 1)
ggsave('../figure/TR.LTR_H3K9me3.line.pdf',width = 3.4,height = 3.3)

LTR_K27 <- data.frame(
  LTR = LTR$X4,
  K27 =H3K27$X4
)
summary(lm(log10(H3K27$X4)~log10(LTR$X4+0.1)))
ggplot(LTR_K27) + geom_abline(slope = 0.83475,intercept = -2.09558,color = Set1[1],size = 1) +
  geom_jitter(aes(x = log10(LTR),y = log10(K27))) + 
  theme_Publication() + 
  labs(x = 'LTR density (log10)',y = 'H3K27me3 intensity (log10)') +
  annotate("text", x = 0.05, y = max(log10(LTR_K9$K9))+1, 
           label = paste("y = 0.83x - 2.09\nR² = 0.21, P < 0.001"),
           hjust = 0, vjust = 1)
ggsave('../figure/TR.LTR_H3K27me3.line.pdf',width = 3.4,height = 3.3)



# Tandem region TR vs. heterochromatin marks ------------------------------
library(reshape2)
library(pheatmap)



# *** 520 INT -------------------------------------------------------------


random <- read_delim('../data/TandemRepeat/random.INT.mafft.matrix',delim = '\t',col_names = T)
random.matrix <- as.matrix((random[,-1]))
pheatmap((random.matrix),show_rownames =F,show_colnames = F,cluster_rows = F,cluster_cols = F,
         filename= '../figure/random.INT.heatmap.pdf',width=3.6,height=3.4)
dev.off()


RT_520 <- read_delim('../data/TandemRepeat/520.INT.mafft.matrix',delim = '\t',col_names = T)
RT_520.matrix <- as.matrix((RT_520[,-1]))
pheatmap((RT_520.matrix),show_rownames =F,show_colnames = F,cluster_rows = F,cluster_cols = F,
         filename= '../figure/RT_520.INT.heatmap.pdf',width=3.6,height=3.4)
dev.off()


random.num <- as.numeric(random.matrix)
RT_520.num <- as.numeric(RT_520.matrix)

x <- max(length(random.num),length(RT_520.num))
length(RT_520.num) <- x
length(random.num) <- x
result <- data.frame(
  RT_520 = RT_520.num,
  random = random.num
) %>% reshape2::melt()

ggplot(result) + stat_ecdf(aes(color = variable,x = value)) + xlim(1,100) + 
  theme_Publication() + scale_colour_Publication('Set1') + 
  labs(x = 'Sequence similarity (%)',y = 'Cumulative fraction',color = 'Type') 
ggsave('../figure/RT_520.INT.ecdf.pdf',width = 2.8,height = 2.8)



# *** 571 INT -------------------------------------------------------------


random <- read_delim('../data/TandemRepeat/random.INT.mafft.matrix',delim = '\t',col_names = T)
random.matrix <- as.matrix((random[,-1]))
pheatmap((random.matrix),show_rownames =F,show_colnames = F,cluster_rows = F,cluster_cols = F,
         filename= '../figure/random.INT.heatmap.pdf',width=3.6,height=3.4)
dev.off()


RT_571 <- read_delim('../data/TandemRepeat/571.INT.mafft.matrix',delim = '\t',col_names = T)
RT_571.matrix <- as.matrix((RT_571[,-1]))
pheatmap((RT_571.matrix),show_rownames =F,show_colnames = F,cluster_rows = F,cluster_cols = F,
         filename= '../figure/RT_571.INT.heatmap.pdf',width=3.6,height=3.4)
dev.off()


random.num <- as.numeric(random.matrix)
RT_571.num <- as.numeric(RT_571.matrix)

x <- max(length(random.num),length(RT_571.num))
length(RT_571.num) <- x
length(random.num) <- x
result <- data.frame(
  RT_571 = RT_571.num,
  random = random.num
) %>% reshape2::melt()

ggplot(result) + stat_ecdf(aes(color = variable,x = value)) + xlim(1,100) + 
  theme_Publication() + scale_colour_Publication('Set1') + 
  labs(x = 'Sequence similarity (%)',y = 'Cumulative fraction',color = 'Type') 
ggsave('../figure/RT_571.INT.ecdf.pdf',width = 2.8,height = 2.8)
# *** 442 INT -------------------------------------------------------------


random <- read_delim('../data/TandemRepeat/random.INT.mafft.matrix',delim = '\t',col_names = T)
random.matrix <- as.matrix((random[,-1]))
pheatmap((random.matrix),show_rownames =F,show_colnames = F,cluster_rows = F,cluster_cols = F,
         filename= '../figure/random.INT.heatmap.pdf',width=3.6,height=3.4)
dev.off()


RT_442 <- read_delim('../data/TandemRepeat/442.INT.mafft.matrix',delim = '\t',col_names = T)
RT_442.matrix <- as.matrix((RT_442[,-1]))
pheatmap((RT_442.matrix),show_rownames =F,show_colnames = F,cluster_rows = F,cluster_cols = F,
         filename= '../figure/RT_442.INT.heatmap.pdf',width=3.6,height=3.4)
dev.off()


random.num <- as.numeric(random.matrix)
RT_442.num <- as.numeric(RT_442.matrix)

x <- max(length(random.num),length(RT_442.num))
length(RT_442.num) <- x
length(random.num) <- x
result <- data.frame(
  RT_442 = RT_442.num,
  random = random.num
) %>% reshape2::melt()

ggplot(result) + stat_ecdf(aes(color = variable,x = value)) + xlim(1,100) + 
  theme_Publication() + scale_colour_Publication('Set1') + 
  labs(x = 'Sequence similarity (%)',y = 'Cumulative fraction',color = 'Type') 
ggsave('../figure/RT_442.INT.ecdf.pdf',width = 2.8,height = 2.8)

# *** 520 RT -------------------------------------------------------------


random <- read_delim('../data/TandemRepeat/random.RT.mafft.matrix',delim = '\t',col_names = T)
random.matrix <- as.matrix((random[,-1]))
pheatmap((random.matrix),show_rownames =F,show_colnames = F,cluster_rows = F,cluster_cols = F,
         filename= '../figure/random.RT.heatmap.pdf',width=3.6,height=3.4)
dev.off()


RT_520 <- read_delim('../data/TandemRepeat/520.RT.mafft.matrix',delim = '\t',col_names = T)
RT_520.matrix <- as.matrix((RT_520[,-1]))
pheatmap((RT_520.matrix),show_rownames =F,show_colnames = F,cluster_rows = F,cluster_cols = F,
         filename= '../figure/RT_520.RT.heatmap.pdf',width=3.6,height=3.4)
dev.off()


random.num <- as.numeric(random.matrix)
RT_520.num <- as.numeric(RT_520.matrix)

x <- max(length(random.num),length(RT_520.num))
length(RT_520.num) <- x
length(random.num) <- x
result <- data.frame(
  RT_520 = RT_520.num,
  random = random.num
) %>% reshape2::melt()

ggplot(result) + stat_ecdf(aes(color = variable,x = value)) + xlim(1,100) + 
  theme_Publication() + scale_colour_Publication('Set1') + 
  labs(x = 'Sequence similarity (%)',y = 'Cumulative fraction',color = 'Type') 
ggsave('../figure/RT_520.RT.ecdf.pdf',width = 2.8,height = 2.8)





# *** 520 RH -------------------------------------------------------------


random <- read_delim('../data/TandemRepeat/random.RH.mafft.matrix',delim = '\t',col_names = T)
random.matrix <- as.matrix((random[,-1]))
pheatmap((random.matrix),show_rownames =F,show_colnames = F,cluster_rows = F,cluster_cols = F,
         filename= '../figure/random.RH.heatmap.pdf',width=3.6,height=3.4)
dev.off()


RT_520 <- read_delim('../data/TandemRepeat/520.RH.mafft.matrix',delim = '\t',col_names = T)
RT_520.matrix <- as.matrix((RT_520[,-1]))
pheatmap((RT_520.matrix),show_rownames =F,show_colnames = F,cluster_rows = F,cluster_cols = F,
         filename= '../figure/RT_520.RH.heatmap.pdf',width=3.6,height=3.4)
dev.off()


random.num <- as.numeric(random.matrix)
RT_520.num <- as.numeric(RT_520.matrix)

x <- max(length(random.num),length(RT_520.num))
length(RT_520.num) <- x
length(random.num) <- x
result <- data.frame(
  RT_520 = RT_520.num,
  random = random.num
) %>% reshape2::melt()

ggplot(result) + stat_ecdf(aes(color = variable,x = value)) + xlim(1,100) + 
  theme_Publication() + scale_colour_Publication('Set1') + 
  labs(x = 'Sequence similarity (%)',y = 'Cumulative fraction',color = 'Type') 
ggsave('../figure/RT_520.RH.ecdf.pdf',width = 2.8,height = 2.8)




# *** 571 RT -------------------------------------------------------------


random <- read_delim('../data/TandemRepeat/random.RT.mafft.matrix',delim = '\t',col_names = T)
random.matrix <- as.matrix((random[,-1]))
pheatmap((random.matrix),show_rownames =F,show_colnames = F,cluster_rows = F,cluster_cols = F,
         filename= '../figure/random.RT.heatmap.pdf',width=3.6,height=3.4)
dev.off()


RT_571 <- read_delim('../data/TandemRepeat/571.RT.mafft.matrix',delim = '\t',col_names = T)
RT_571.matrix <- as.matrix((RT_571[,-1]))
pheatmap((RT_571.matrix),show_rownames =F,show_colnames = F,cluster_rows = F,cluster_cols = F,
         filename= '../figure/RT_571.RT.heatmap.pdf',width=3.6,height=3.4)
dev.off()


random.num <- as.numeric(random.matrix)
RT_571.num <- as.numeric(RT_571.matrix)

x <- max(length(random.num),length(RT_571.num))
length(RT_571.num) <- x
length(random.num) <- x
result <- data.frame(
  RT_571 = RT_571.num,
  random = random.num
) %>% reshape2::melt()

ggplot(result) + stat_ecdf(aes(color = variable,x = value)) + xlim(1,100) + 
  theme_Publication() + scale_colour_Publication('Set1') + 
  labs(x = 'Sequence similarity (%)',y = 'Cumulative fraction',color = 'Type') 
ggsave('../figure/RT_571.RT.ecdf.pdf',width = 2.8,height = 2.8)





# *** 571 RH -------------------------------------------------------------


random <- read_delim('../data/TandemRepeat/random.RH.mafft.matrix',delim = '\t',col_names = T)
random.matrix <- as.matrix((random[,-1]))
pheatmap((random.matrix),show_rownames =F,show_colnames = F,cluster_rows = F,cluster_cols = F,
         filename= '../figure/random.RH.heatmap.pdf',width=3.6,height=3.4)
dev.off()


RT_571 <- read_delim('../data/TandemRepeat/571.RH.mafft.matrix',delim = '\t',col_names = T)
RT_571.matrix <- as.matrix((RT_571[,-1]))
pheatmap((RT_571.matrix),show_rownames =F,show_colnames = F,cluster_rows = F,cluster_cols = F,
         filename= '../figure/RT_571.RH.heatmap.pdf',width=3.6,height=3.4)
dev.off()


random.num <- as.numeric(random.matrix)
RT_571.num <- as.numeric(RT_571.matrix)

x <- max(length(random.num),length(RT_571.num))
length(RT_571.num) <- x
length(random.num) <- x
result <- data.frame(
  RT_571 = RT_571.num,
  random = random.num
) %>% reshape2::melt()

ggplot(result) + stat_ecdf(aes(color = variable,x = value)) + xlim(1,100) + 
  theme_Publication() + scale_colour_Publication('Set1') + 
  labs(x = 'Sequence similarity (%)',y = 'Cumulative fraction',color = 'Type') 
ggsave('../figure/RT_571.RH.ecdf.pdf',width = 2.8,height = 2.8)





# *** 442 RT -------------------------------------------------------------


random <- read_delim('../data/TandemRepeat/random.RT.mafft.matrix',delim = '\t',col_names = T)
random.matrix <- as.matrix((random[,-1]))
pheatmap((random.matrix),show_rownames =F,show_colnames = F,cluster_rows = F,cluster_cols = F,
         filename= '../figure/random.RT.heatmap.pdf',width=3.6,height=3.4)
dev.off()


RT_442 <- read_delim('../data/TandemRepeat/442.RT.mafft.matrix',delim = '\t',col_names = T)
RT_442.matrix <- as.matrix((RT_442[,-1]))
pheatmap((RT_442.matrix),show_rownames =F,show_colnames = F,cluster_rows = F,cluster_cols = F,
         filename= '../figure/RT_442.RT.heatmap.pdf',width=3.6,height=3.4)
dev.off()


random.num <- as.numeric(random.matrix)
RT_442.num <- as.numeric(RT_442.matrix)

x <- max(length(random.num),length(RT_442.num))
length(RT_442.num) <- x
length(random.num) <- x
result <- data.frame(
  RT_442 = RT_442.num,
  random = random.num
) %>% reshape2::melt()

ggplot(result) + stat_ecdf(aes(color = variable,x = value)) + xlim(1,100) + 
  theme_Publication() + scale_colour_Publication('Set1') + 
  labs(x = 'Sequence similarity (%)',y = 'Cumulative fraction',color = 'Type') 
ggsave('../figure/RT_442.RT.ecdf.pdf',width = 2.8,height = 2.8)





# *** 442 RH -------------------------------------------------------------


random <- read_delim('../data/TandemRepeat/random.RH.mafft.matrix',delim = '\t',col_names = T)
random.matrix <- as.matrix((random[,-1]))
pheatmap((random.matrix),show_rownames =F,show_colnames = F,cluster_rows = F,cluster_cols = F,
         filename= '../figure/random.RH.heatmap.pdf',width=3.6,height=3.4)
dev.off()


RT_442 <- read_delim('../data/TandemRepeat/442.RH.mafft.matrix',delim = '\t',col_names = T)
RT_442.matrix <- as.matrix((RT_442[,-1]))
pheatmap((RT_442.matrix),show_rownames =F,show_colnames = F,cluster_rows = F,cluster_cols = F,
         filename= '../figure/RT_442.RH.heatmap.pdf',width=3.6,height=3.4)
dev.off()


random.num <- as.numeric(random.matrix)
RT_442.num <- as.numeric(RT_442.matrix)

x <- max(length(random.num),length(RT_442.num))
length(RT_442.num) <- x
length(random.num) <- x
result <- data.frame(
  RT_442 = RT_442.num,
  random = random.num
) %>% reshape2::melt()

ggplot(result) + stat_ecdf(aes(color = variable,x = value)) + xlim(1,100) + 
  theme_Publication() + scale_colour_Publication('Set1') + 
  labs(x = 'Sequence similarity (%)',y = 'Cumulative fraction',color = 'Type') 
ggsave('../figure/RT_520.RH.ecdf.pdf',width = 2.8,height = 2.8)




# RNAfold -----------------------------------------------------------------


seq442 <- read_delim('../data/TandemRepeat/RNA_fold/442.w300_step1.value',delim = '\t',col_names = F) %>% 
  mutate(Type = 'pSc119.2')
seq380 <- read_delim('../data/TandemRepeat/RNA_fold/380.w300_step1.value',delim = '\t',col_names = F)%>% 
  mutate(Type = 'pSc200')
seq571 <- read_delim('../data/TandemRepeat/RNA_fold/571.w300_step1.value',delim = '\t',col_names = F)%>% 
  mutate(Type = 'pSc250')
seqRandom <- read_delim('../data/TandemRepeat/RNA_fold/random.seed2.value',delim = '\t',col_names = F)%>% 
  mutate(Type = 'Random')
result <- rbind(seq380,seq442,seq571,seqRandom)



ggplot(result,aes(x = Type, y = X1,fill = Type)) + geom_boxplot(size = 0.2,outlier.alpha = 0.1) + 
  scale_fill_Publication() + theme_Publication() + 
  labs(x = 'Type', y = 'Gibbs free energy (kcal/mol)') + 
  theme(
    legend.position = 'None'
  ) + 
  geom_signif(comparisons = list(c("pSc119.2", "Random")), y_position = -20,
              map_signif_level=TRUE)+ 
  geom_signif(comparisons = list(c("pSc200", "Random")), y_position = -13,
              map_signif_level=TRUE)+ 
  geom_signif(comparisons = list(c("pSc250", "Random")), y_position = -7,
              map_signif_level=TRUE)



ggsave('../figure/RNAfold.compare.boxplot.pdf',width = 3.65,height = 3.23)
  