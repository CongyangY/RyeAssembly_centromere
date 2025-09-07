library(tidyverse)
library(Biostrings)


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

# 1. cal PCA -----------------------------------------------------------------

aligned.seq <- readDNAStringSet('project/WN_YCY/result/anno/TE/manual_repeat.520.fas')
aligned.seq <- as.data.frame(aligned.seq)
seq_matrix <-sapply(t(tmp),function(x) {
  sapply(strsplit(x, ""), function(y) {
    unlist(lapply(y, function(nucleotide) {
      switch(nucleotide,
             "A" = 1,
             "C" = 2,
             "T" = 3,
             "G" = 4,
             "-" = 0)
    }))
  })
})
colnames(seq_matrix) <- rownames(aligned.seq)
seq_matrix <- t(seq_matrix)  # 转置矩阵以匹配PCA函数的需求

pca_result <- prcomp(seq_matrix,scale. = T)
pca_data <- data.frame(PC1 = pca_result$x[,1], PC2 = pca_result$x[,2], RowNames = rownames(pca_result$x))

ggplot(pca_data, aes(x = PC1, y = PC2)) +
  geom_point()+
  xlab("PC1") +
  ylab("PC2") +
  ggtitle("PCA of Aligned Sequences")




var_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2)
pca_screen <- data.frame(PC = 1:length(var_explained), Variance = var_explained)
ggplot(pca_screen[1:10,], aes(x = PC, y = Variance)) +
  geom_point() + # 添加点
  geom_line() + # 添加线
  scale_x_continuous(breaks = 1:length(var_explained)) + # 调整X轴刻度
  xlab("Principal Component") + # X轴标签
  ylab("Variance Explained") + # Y轴标签
  ggtitle("Scree Plot") 

write_delim(pca_screen,'../data/repeat_consensus/manual_repeat.520.pca.screen',delim = '\t')





scores <- as.data.frame(pca_result$x[, 1:2])

# 应用K-均值聚类
set.seed(123) # 确保结果可重复
clusters <- kmeans(scores, centers = 3) 
scores$cluster <- as.factor(clusters$cluster)
ggplot(scores, aes(x = PC1, y = PC2, color = cluster)) +
  geom_point() +
  theme_minimal() +
  ggtitle("Classification based on PC1 and PC2") +
  xlab("PC1") +
  ylab("PC2")

scores$Rownames = rownames(pca_result$x)
write_delim(scores,'project/WN_YCY/result/anno/TE/manual_repeat.520_1.cluster.pca',delim = '\t')

write_delim(filter(scores,cluster == 1),
            'project/WN_YCY/result/anno/TE/manual_repeat.520_1.list',delim = '\t')

write_delim(filter(scores,cluster == 2),
            'project/WN_YCY/result/anno/TE/manual_repeat.520_2.list',delim = '\t')
write_delim(filter(scores,cluster == 3),
            'project/WN_YCY/result/anno/TE/manual_repeat.520_3.list',delim = '\t')


# 2.  plot  PCA and screen -----------------------------------------------------------------

screen.442 <- read_delim('../data/repeat_consensus/manual_repeat.442.pca.screen',delim = '\t')
ggplot(screen.442[1:10,], aes(x = PC, y = Variance)) +
  geom_point() + # 添加点
  geom_line() + # 添加线
  scale_x_continuous(breaks = 1:10) + # 调整X轴刻度
  xlab("Principal Component") + # X轴标签
  ylab("Variance Explained") + # Y轴标签
  ggtitle("Scree Plot") + theme_Publication()
ggsave('../figure/repeat_consensus/tr.442.pca.screen.pdf',width = 6,height = 4)

pca_data.442 <- read_delim('../data/repeat_consensus/manual_repeat.442.pca.cluster',delim = '\t')
pca_data.442$cluster <- as.factor(pca_data.442$cluster)
ggplot(pca_data.442, aes(x = PC1, y = PC2, color = cluster)) +
  geom_point() +
  theme_minimal() +
  ggtitle("Cluster based on PC1 and PC2") +
  xlab("PC1") +
  ylab("PC2") + theme_Publication() + scale_colour_Publication('Set2')
ggsave('../figure/repeat_consensus/tr.442.pca.cluster.pdf',width =8 ,height = 6)






screen.520 <- read_delim('../data/repeat_consensus/manual_repeat.520.pca.screen',delim = '\t')
ggplot(screen.520[1:10,], aes(x = PC, y = Variance)) +
  geom_point() + # 添加点
  geom_line() + # 添加线
  scale_x_continuous(breaks = 1:10) + # 调整X轴刻度
  xlab("Principal Component") + # X轴标签
  ylab("Variance Explained") + # Y轴标签
  ggtitle("Scree Plot") + theme_Publication()
ggsave('../figure/repeat_consensus/tr.520.pca.screen.pdf',width = 6,height = 4)

pca_data.520<- read_delim('../data/repeat_consensus/manual_repeat.520.pca.cluster',delim = '\t')
pca_data.520$cluster <- as.factor(pca_data.520$cluster)
ggplot(pca_data.520, aes(x = PC1, y = PC2, color = cluster)) +
  geom_point() +
  theme_minimal() +
  ggtitle("Cluster based on PC1 and PC2") +
  xlab("PC1") +
  ylab("PC2") + theme_Publication() + scale_colour_Publication('Set2')
ggsave('../figure/repeat_consensus/tr.520.pca.cluster.pdf',width =8 ,height = 6)


