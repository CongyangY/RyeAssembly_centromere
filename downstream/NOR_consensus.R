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


# length distribution -----------------------------------------------------

len.dis <- read_delim('../data/NOR/rye_45S.unit.bed',delim = '\t',col_names = F) %>% 
  mutate(len = X3-X2,color = ifelse(len == 5809,'max','other'))

ggplot(len.dis) + geom_bar(aes(x = len,color = color,fill = color)) + 
  xlim(c(5700,5900))

# # 1. cal PCA -----------------------------------------------------------------
# 
# aligned.seq <- readDNAStringSet('../data/NOR/rye_45S.unit.mafft.fa')
# aligned.seq <- as.data.frame(aligned.seq)
# seq_matrix <-sapply(t(aligned.seq),function(x) {
#   sapply(strsplit(x, ""), function(y) {
#     unlist(lapply(y, function(nucleotide) {
#       switch(nucleotide,
#              "A" = 1,
#              "C" = 2,
#              "T" = 3,
#              "G" = 4,
#              "-" = 0)
#     }))
#   })
# })
# colnames(seq_matrix) <- rownames(aligned.seq)
# seq_matrix <- t(seq_matrix)  # 转置矩阵以匹配PCA函数的需求
# 
# pca_result <- prcomp(seq_matrix,scale. = T)
# pca_data <- data.frame(PC1 = pca_result$x[,1], PC2 = pca_result$x[,2], RowNames = rownames(pca_result$x))
# 
# ggplot(pca_data, aes(x = PC1, y = PC2)) +
#   geom_point()+
#   xlab("PC1") +
#   ylab("PC2") +
#   ggtitle("PCA of Aligned Sequences")
# 
# 
# 
# 
# var_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2)
# pca_screen <- data.frame(PC = 1:length(var_explained), Variance = var_explained)
# ggplot(pca_screen[1:10,], aes(x = PC, y = Variance)) +
#   geom_point() + # 添加点
#   geom_line() + # 添加线
#   scale_x_continuous(breaks = 1:length(var_explained)) + # 调整X轴刻度
#   xlab("Principal Component") + # X轴标签
#   ylab("Variance Explained") + # Y轴标签
#   ggtitle("Scree Plot") 
# 
# 
# 
# scores <- as.data.frame(pca_result$x[, 1:2])
# 
# # 应用K-均值聚类
# set.seed(123) # 确保结果可重复
# clusters <- kmeans(scores, centers = 2) 
# scores$cluster <- as.factor(clusters$cluster)
# ggplot(scores, aes(x = PC1, y = PC2, color = cluster)) +
#   geom_point() +
#   theme_minimal() +
#   ggtitle("Classification based on PC1 and PC2") +
#   xlab("PC1") +
#   ylab("PC2")
# 
# scores$Rownames = rownames(pca_result$x)
# write_delim(scores,'../data/NOR/rye_45S_1.unit.fas',delim = '\t')
# 
# write_delim(filter(scores,cluster == 1),
#             '../data/NOR/rye_45S_1.list',delim = '\t')
# 
# write_delim(filter(scores,cluster == 2),
#             'project/WN_YCY/result/anno/TE/manual_repeat.520_2.list',delim = '\t')
# write_delim(filter(scores,cluster == 3),
#             'project/WN_YCY/result/anno/TE/manual_repeat.520_3.list',delim = '\t')
# 
# 
# 

# pos variation -----------------------------------------------------------

pos_cov.NOR <- read_delim('../data/NOR/rye_45S.unit.5809.pos_cov',delim = '\t',col_names = F)


ggplot(pos_cov.NOR) + geom_bar(aes(x = X1,y=X2),stat = 'identity')  + ylim(c(0,0.1)) + 
  theme_Publication() + labs(x = 'Position (bp)', y = 'Variant frequency')
ggsave('../figure/NOR/NOR_45s_rDNA.pos_cov.barplot.pdf',width = 6, height = 4)



# methylation level -------------------------------------------------------------

aver <- function(data){mean(data,na.rm = T,trim = 0.0025)}

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


CG.NOR <- plot_data('../data/NOR/CG.NOR.tandemRepeat.matrix.gz','45s')
CHG.NOR <- plot_data('../data/NOR/CHG.NOR.tandemRepeat.matrix.gz','45s')
CHH.NOR <- plot_data('../data/NOR/CHH.NOR.tandemRepeat.matrix.gz','45s')

ggplot(data=CG.NOR,aes(x=pos)) + 
  geom_line(aes(group = variable,y = value,color = variable),size = 0.8) + 
  labs(x = '',
       y = 'CG methylaton level\n',
       color = '')   +
  scale_x_continuous(
    breaks = c(1,30,90,120),
    labels = c('-3K','left','right','3K'),
    expand = c(0.005,0.005)) + 
  theme_bw() + theme_Publication() +
  scale_colour_Publication('Accent') + theme(
    legend.position = 'None'
  )



# methylation position ----------------------------------------------------
CG <- read_delim('../data/NOR/rye_45S.unit.CG.pos.bed',delim = '\t',col_names = F)
CHG <- read_delim('../data/NOR/rye_45S.unit.CHG.pos.bed',delim = '\t',col_names = F)
CHH <- read_delim('../data/NOR/rye_45S.unit.CHH.pos.bed',delim = '\t',col_names = F)


CG.stat <- as.data.frame(table(cut(CHG$X2,breaks = c(100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500,1600,1700,1800,
                                        1900,2000,2100,2200,2300,2400,2500,2600,2700,2800,2900,3000,3100,3200,3300,3400,
                                        3500,3600,3700,3800,3900,4000,4100,4200,4300,4400,4500,4600,4700,4800,4900,5000,
                                        5100,5200,5300,5400,5500,5600,5700,5809),
                        labels = c('100-200', '200-300', '300-400', '400-500', '500-600',
                                   '600-700', '700-800', '800-900', '900-1000', '1000-1100',
                                   '1100-1200', '1200-1300', '1300-1400', '1400-1500',
                                   '1500-1600', '1600-1700', '1700-1800', '1800-1900',
                                   '1900-2000', '2000-2100', '2100-2200', '2200-2300',
                                   '2300-2400', '2400-2500', '2500-2600', '2600-2700',
                                   '2700-2800', '2800-2900', '2900-3000', '3000-3100',
                                   '3100-3200', '3200-3300', '3300-3400', '3400-3500',
                                   '3500-3600', '3600-3700', '3700-3800', '3800-3900',
                                   '3900-4000', '4000-4100', '4100-4200', '4200-4300',
                                   '4300-4400', '4400-4500', '4500-4600', '4600-4700',
                                   '4700-4800', '4800-4900', '4900-5000', '5000-5100',
                                   '5100-5200', '5200-5300', '5300-5400', '5400-5500',
                                   '5500-5600', '5600-5700', '>5700'))))

ggplot(CG) + geom_bar(aes(x = X2),bins)



start <- 100
end <- 5800
step <- 10
breaks <- seq(from = start, to = end + step, by = step)
labels <- sapply(1:(length(breaks)-1), function(i) paste(breaks[i], breaks[i+1]-1, sep="-"))
labels[length(labels)] <- paste(">5800")
CG.stat <- as.data.frame(table(cut(CHH$X2, breaks = breaks, labels = labels, right = FALSE)))

ggplot(CG.stat) + geom_bar(aes(x = Var1, y = Freq),stat = 'identity') + 
  labs(x = '')  + theme_Publication()+ theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )
