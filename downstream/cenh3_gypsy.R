library(tidyverse)

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


# 1. Gypsy vs. CENH3 ------------------------------------------------------


cenh3 <- read_delim('../data/centromere/WN.pericent40M.w1M.CENH3_mean.bg',
                    delim = '\t',col_names = F) %>% mutate(
                      Type = 'CENH3'
                    )
Gypsy <- read_delim('../data/centromere/Gypsy.intact.TEcontent.bg',
                    delim = '\t',col_names = F) %>% mutate(
                      Type = 'Gypsy'
                    )
InsertTime <- read_delim('../data/centromere/Gypsy.intact.insertTime.mean.bg',
                    delim = '\t',col_names = F) %>% mutate(
                      Type = 'InsertTime'
                    )

result <- data.frame(
  CENH3 = cenh3$X4,
  Gypsy = Gypsy$X4/1000000,
  InsertTime = InsertTime$X4
)
result <- reshape(result,id = 'Type')
ggplot(result) + 
  geom_jitter(aes(x = CENH3, y = Gypsy,size = InsertTime,color = InsertTime),
             alpha = 0.15) + 
  scale_color_gradient(low = "blue", high = "red") + 
  theme_Publication() + theme(
    legend.direction = 'vertical',
    legend.position = 'right'
  )
ggsave('../figure/Gypsy_CENH3.dotplot.pdf',width = 4.35,height = 3.5)


library(reshape2)
cor_matrix <- cor(result[c("CENH3", "Gypsy", "InsertTime")])

# 将矩阵转换为长格式
cor_melted <- melt(cor_matrix)

# 创建热图
ggplot(cor_melted, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-1,1)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "", y = "", fill = "Correlation") +
  geom_text(aes(label = sprintf("%.2f", value)), color = "black", size = 4)

ggsave('../figure/Gypsy_CENH3.corplot.pdf',width = 4.35,height = 3.15)




# Ka/Ks -------------------------------------------------------------------

chr1R.cent <- read_delim('../data/centromere/cent_evol/Gypsy.cent.RT.chr1R.ka_ks',
                         delim = ',',skip = 1,col_names = c('X1','X2','Ratio')) %>% 
  mutate(Ratio = as.numeric(Ratio),Type = 'Cent')
chr1R.shuffle <- read_delim('../data/centromere/cent_evol/Gypsy.young.no_cent.RT.chr1R.ka_ks',
                         delim = ',',skip = 1,col_names = c('X1','X2','Ratio')) %>% 
  mutate(Ratio = as.numeric(Ratio),Type = 'Shuffle')


chr1R <- rbind(chr1R.cent,chr1R.shuffle)
ggplot(chr1R) + geom_boxplot(aes(x = Type, y = Ratio,fill = Type))

ggplot(chr1R) + geom_density(aes(x = Ratio, color = Type))



cent <- read_delim('../data/centromere/cent_evol/Gypsy.cent.RT.ka_ks',
                   delim = ',',skip = 1,col_names = c('X1','X2','Ratio')) %>% 
  mutate(Ratio = as.numeric(Ratio),Type = 'Cent')
shuffle <- read_delim('../data/centromere/cent_evol/Gypsy.young.no_cent.RT.ka_ks',
                   delim = ',',skip = 1,col_names = c('X1','X2','Ratio')) %>% 
  mutate(Ratio = as.numeric(Ratio),Type = 'Shuffle')


library(ggsignif)
result <- rbind(cent,shuffle)
ggplot(result,aes(x = Type, y = Ratio)) + geom_violin() + 
  geom_signif(comparisons=list(c("Cent", "Shuffle")),map_signif_level = T) 
  
ggplot(result,aes(color = Type, x = Ratio)) + stat_ecdf() + 
  theme_Publication() + scale_colour_Publication('Set1')
ggsave('../figure/Ka_Ks.Gypsy.ecdf.pdf',width = 3.6,height = 3.2)  



# Gypsy CRM ---------------------------------------------------------------
library('pheatmap')

gypsy.rpkm <- read_delim('../data/centromere/cent_evol/rye.ChIP.Gypsy.RPKM',
                         delim = '\t',col_names = T)
gypsy.rpkm.matrix <- as.matrix(gypsy.rpkm[,-1])
rownames(gypsy.rpkm.matrix) <- as.character(gypsy.rpkm$Gypsy)
pheatmap(gypsy.rpkm.matrix,cluster_col = T,cluster_rows = F,show_colnames = T,
         show_rownames = T,width=3,height=5,angle_col = "45",display_numbers = TRUE)


# CRM ---------------------------------------------------------------------

WW.ChIP <- read_delim('../data/centromere/cent_evol/WW-ChIP.CRM.w100_s50.bg',delim = '\t',col_names = F) %>% 
  mutate(Type = 'WW')
WN.ChIP <- read_delim('../data/centromere/cent_evol/WN-ChIP.CRM.w100_s50.bg',delim = '\t',col_names = F) %>% 
  mutate(Type = 'WN')
IMP.ChIP <- read_delim('../data/centromere/cent_evol/IMP-ChIP.CRM.w100_s50.bg',delim = '\t',col_names = F) %>% 
  mutate(Type = 'IMP')
SD.ChIP <- read_delim('../data/centromere/cent_evol/SD-ChIP.CRM.w100_s50.bg',delim = '\t',col_names = F) %>% 
  mutate(Type = 'SD')
LD.ChIP <- read_delim('../data/centromere/cent_evol/LD-ChIP.CRM.w100_s50.bg',delim = '\t',col_names = F) %>% 
  mutate(Type = 'LD')
Lo7.ChIP <- read_delim('../data/centromere/cent_evol/Lo7-ChIP.CRM.w100_s50.bg',delim = '\t',col_names = F) %>% 
  mutate(Type = 'Lo7')
tmp <- rbind(WW.ChIP,WN.ChIP,IMP.ChIP,SD.ChIP,LD.ChIP,Lo7.ChIP)
ggplot(tmp) + geom_line(aes(x = X2,y = X4,color = Type))


WW <- read_delim('../data/centromere/cent_evol/WW-CENH3.wheat_CRM.w40.bg',delim = '\t',col_names = F) %>% 
  mutate(Type = 'WW',X4 = ifelse(X4<0,0,X4))
WN <- read_delim('../data/centromere/cent_evol/WN-CENH3.wheat_CRM.w40.bg',delim = '\t',col_names = F) %>% 
  mutate(Type = 'WN',X4 = ifelse(X4<0,0,X4))
IMP <- read_delim('../data/centromere/cent_evol/IMP-ChIP.wheat_CRM.w40.bg',delim = '\t',col_names = F) %>% 
  mutate(Type = 'IMP',X4 = ifelse(X4<0,0,X4))
SD <- read_delim('../data/centromere/cent_evol/SD-ChIP.wheat_CRM.w40.bg',delim = '\t',col_names = F) %>% 
  mutate(Type = 'SD',X4 = ifelse(X4<0,0,X4))
LD <- read_delim('../data/centromere/cent_evol/LD-CENH3.wheat_CRM.w40.bg',delim = '\t',col_names = F) %>% 
  mutate(Type = 'LD',X4 = ifelse(X4<0,0,X4))
Lo7 <- read_delim('../data/centromere/cent_evol/Lo7-ChIP.wheat_CRM.w40.bg',delim = '\t',col_names = F) %>% 
  mutate(Type = 'Lo7',X4 = ifelse(X4<0,0,X4))
wheat_CRM <- rbind(WW,WN,IMP,SD,LD,Lo7)
ggplot(tmp) + geom_line(aes(x = X2,y = X4,color = Type)) + 
  scale_colour_Publication() + theme_Publication() + 
  labs(x = 'CRM Postion (bp)', y = 'Relative CENH3 signal')
ggsave('../figure/CRM.CENH3_binding.metaplot.pdf',width = 4.2,height = 3)

wheat_CRM.CV <- data.frame(
  WN = WN$X4,
  IMP = IMP$X4,
  Lo7 = Lo7$X4,
  WW = WW$X4,
  LD = LD$X4,
  SD = SD$X4
)
CV_cal <- function(data){
  return(sd(data[1:6])/mean(data[1:6]))
}

p_cal <- function(data){
  return(t.test(data[1:3],data[5:6])$p.value)
}
fc_cal <- function(data){
  return(sum(data[1:3])/sum(data[5:6]))
}
wheat_CRM.CV$p = apply(wheat_CRM.CV,1,p_cal)
wheat_CRM.CV$fc = apply(wheat_CRM.CV,1,fc_cal)
wheat_CRM.CV$p <- -log10(wheat_CRM.CV$p)
wheat_CRM.CV$p_scale = -log10(0.05)
wheat_CRM.CV$fc_scale = 2
p_value_colors <- colorRampPalette(rev(c("#67001f", "#b2182b", "#d6604d", "#f4a582", "#fddbc7", "#f7f7f7")))(100)
pheatmap(wheat_CRM.CV[,c(7,9)],cluster_rows = F,
         width=2.3,height=3,show_rownames = F, show_colnames = F,
         filename = '../figure/wheat_CRM.p_value.heatmap.pdf',color = p_value_colors)
fc_colors <- colorRampPalette(c("#00441b", "#1b7837", "#5aae61", "#a6dba0", "#d9f0d3", "#f7f7f7", "#e7d4e8", "#c2a5cf", "#9970ab", "#762a83", "#40004b"))(100)
pheatmap(wheat_CRM.CV[,c(8,10)],cluster_rows = F,
         width=2.3,height=3,show_rownames = F, show_colnames = F,
         filename = '../figure/wheat_CRM.FC.heatmap.pdf',color = fc_colors)



# LTR model ---------------------------------------------------------------

library(ggplot2)
library(ggrepel)

# 创建数据框
domains <- data.frame(
  start = c(2951, 3981, 4786, 5498, 6149, 7160),
  end = c(3331, 4205, 5304, 6028, 7072, 7321),
  name = c("GAG", "PROT", "RT", "RH", "INT", "CHDCR"),
  y = rep(1, 6)
)

# 创建图形
ggplot() +
  # 绘制全长转座子
  geom_segment(aes(x = 0, xend = 8161, y = 1, yend = 1), 
               arrow = arrow(ends = "both", length = unit(0.2, "cm")), 
               size = 1) +
  # 绘制各个domain
  geom_rect(data = domains, 
            aes(xmin = start, xmax = end, ymin = 0.9, ymax = 1.1, fill = name),
            color = "black") +
  # 添加domain标签
  geom_text_repel(data = domains, 
                  aes(x = (start + end) / 2, y = 1.1, label = name),
                  direction = "y", nudge_y = 0.2) +
  # 设置x轴刻度和标签
  scale_x_continuous(breaks = c(0, domains$start, domains$end, 8161),
                     labels = c(0, domains$start, domains$end, 8161)) +
  # 设置主题
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        panel.grid = element_blank()) +
  # 设置标题和坐标轴标签
  labs(title = "Transposon Model (8161bp)",
       x = "Position (bp)",
       fill = "Domain") +
  # 调整填充颜色
  scale_fill_brewer(palette = "Set3")

ggsave("../figure/wheat_CRM.transposon_model.pdf", width = 15, height = 3)
