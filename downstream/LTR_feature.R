library(tidyverse)
library(reshape2)
library(RColorBrewer)
library(pheatmap)
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



# LTR insertion time ------------------------------------------------------


LTR <- read_delim('../data/TE/WNYCY.fa.mod.LTR.intact.LTR_identity.bed',delim = '\t',col_names = F) %>% 
  mutate(insertTime = (1-X5)/2.6e-8)

ggplot(LTR) + geom_density(aes(x = insertTime)) + theme_Publication() + 
  labs(x = 'LTR insertion time (year)',y = 'Density') + 
  geom_vline(xintercept = 3.5e5,color = 'gray',linetype = 2) + 
  geom_vline(xintercept = 1.2e6,color = 'gray',linetype = 2)

ggsave('../figure/circos/LTR_insertionTime.density.pdf',width = 6,height = 4)
write_delim(filter(LTR,insertTime<3.5e5)[,1:6],'../data/TE/LTR_young.bed',delim = '\t',col_names = F)
write_delim(filter(LTR,insertTime>=3.5e5,insertTime<1.2e6)[,1:6],'../data/TE/LTR_mid.bed',delim = '\t',col_names = F)
write_delim(filter(LTR,insertTime>=1.2e6)[,1:6],'../data/TE/LTR_old.bed',delim = '\t',col_names = F)



# LTR length distribution -------------------------------------------------

L_SINE <- read_delim('../data/TE/SINE_LINE.length.chr.list',delim = ' ',col_names = F) %>% 
                       mutate(Type = 'Others',Chr = X1,Length = X2)
transposon <- read_delim('../data/TE/transposon.length.chr.list',delim = ' ',col_names = F) %>% 
  mutate(Type = 'DNA transposon',Chr = X1,Length = X2)
retrotransposon <- read_delim('../data/TE/retrotransposon.length.chr.list',delim = ' ',col_names = F) %>% 
  mutate(Type = 'Retrotransposon',Chr = X1,Length = X2)


TE.dis <- rbind(L_SINE,transposon,retrotransposon)
TE.dis$Type <- factor(TE.dis$Type,levels = c('Others','DNA transposon','Retrotransposon'))
ggplot(TE.dis) + geom_bar(aes(x = Chr,y = Length,color = Type, fill = Type),
                          stat = 'identity',position = 'fill') +
  theme_Publication() + labs(x = 'Chromosome', y = 'Proportion') + 
  scale_fill_Publication('Set1') + scale_colour_Publication('Set1')

ggsave('../figure/circos/LTR_len.dis.bar.svg',width = 6,height = 4)
# cent LTR identity distribution -------------------------------------------------

cent.LTR. <- read_delim('../data/CENH3/LTR.intact.cent.bed',delim = '\t',col_names = F)
ggplot(cent.LTR) + geom_density(aes(x = X5))


aver <- function(data){mean(data,na.rm = T,trim = 0.0025)}
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


LTR.old <- plot_data('../data/CENH3/LTR.intact.cent.low.CENH3.matrix.gz','Old')
LTR.mid<- plot_data('../data/CENH3/LTR.intact.cent.mid.CENH3.matrix.gz','Mid')
LTR.high <- plot_data('../data/CENH3/LTR.intact.cent.high.CENH3.matrix.gz','Young')

LTR.age <- rbind(LTR.old,LTR.mid,LTR.high)

ggplot(data=LTR.age,aes(x=pos)) + 
  geom_line(aes(group = variable,y = value,color = variable),size = 1) + 
  labs(x = '',
       y = 'Mean relative signal\n',
       color = '') +
  scale_x_continuous(
    breaks = c(1,90,150,240),
    labels = c('-6K','TSS','TES','6K'),
    expand = c(0.005,0.005)) + 
  theme_Publication() + scale_colour_Publication('Set2')

ggsave('../figure/CENH3/CENH3.LTR.meta.pdf',width = 5,height = 4)



# LTR vs. nuclear  --------------------------------------------------------

nucl <- read_delim('../data/CENH3/WNYCY.cent.1M.cenh3_nuclear.count',
                   delim = '\t',col_names = F)
LTR <- read_delim('../data/CENH3/WNYCY.cent.1M.flLTR.identity',
                  delim = '\t',col_names = F)
result <- data.frame(
  nucl = nucl$X4,
  LTR = LTR$X4)


ggplot(result) + geom_point(aes(x = nucl,y=LTR))



# LTR family --------------------------------------------------------------

LTR.fam <- read_delim('../data/TE/WNYCY.cent.LTR.family.chr.list',
                      delim = '\t',col_names = F) %>% 
  group_by(X2) %>% 
  mutate(X3 = ifelse(X1<100,'Others',X3))

others <- filter(LTR.fam,X3 == 'Others') %>% 
  group_by(X2) %>% summarise(X1 = sum(X1)) %>% 
  mutate(
    X3 = 'Others',
    X4 = 'Others'
  ) %>% select(X1,X2,X3,X4)
LTR.fam <- rbind(filter(LTR.fam,X3!='Others'),others)  



ggplot(LTR.fam) +geom_bar(aes(x = X2,y=X1,fill = X3,color = X3),stat = 'identity') + 
  theme_Publication() + scale_fill_Publication('Paired') + 
  scale_colour_Publication('Paired') + 
  labs(x = 'Chromosome', y = 'Count',fill = 'Type', color = 'Type') + 
  theme(
    legend.position = 'right',
    legend.direction = 'vertical',
    axis.text.x = element_text(angle = 60,hjust = 1)
  )
ggsave('../figure/CENH3/CENH3.LTR.count.bar.pdf',
       width = 5,height = 4)


# LTR insertion time heatmap ----------------------------------------------


LTR <- read_delim('../data/TE/WNYCY.fa.mod.LTR.intact.LTR_identity.cent.bed',
                  delim = '\t',col_names = F) %>% 
  mutate(
    X2 = ifelse(X1 == 'chr1R',X2-429700000,X2),
    X3 = ifelse(X1 == 'chr1R',X3-429700000,X3),
    
    X2 = ifelse(X1 == 'chr2R',X2-474000000,X2),
    X3 = ifelse(X1 == 'chr2R',X3-474000000,X3),
    
    X2 = ifelse(X1 == 'chr3R',X2-537000000,X2),
    X3 = ifelse(X1 == 'chr3R',X3-537000000,X3),
    
    X2 = ifelse(X1 == 'chr4R',X2-387000000,X2),
    X3 = ifelse(X1 == 'chr4R',X3-387000000,X3),
    
    X2 = ifelse(X1 == 'chr5R',X2-460000000,X2),
    X3 = ifelse(X1 == 'chr5R',X3-460000000,X3),
    
    X2 = ifelse(X1 == 'chr6R',X2-366000000,X2),
    X3 = ifelse(X1 == 'chr6R',X3-366000000,X3),
    
    X2 = ifelse(X1 == 'chr7R',X2-471000000,X2),
    X3 = ifelse(X1 == 'chr7R',X3-471000000,X3),
    
    SNP = paste(X1,X2,sep = '_'),
    Chromosome = str_replace(X1,'chr',''),
    Position = (X2+X3)/2,
    Trait = (1-X5)/2.6e-8
    ) %>% select(
      SNP,Chromosome,Position,Trait
    )
CMplot(LTR,type="p",plot.type="d",bin.size=1e6,chr.den.col=c("darkgreen", "yellow", "red"),file="jpg",
       dpi=300,main="illumilla_60K",file.output=TRUE,verbose=TRUE,width=9,height=6)


pheatmap()

LTR$Chromosome <- as.factor(LTR$Chromosome)
LTR$Position <- as.numeric(LTR$Position)
LTR$Trait <- as.numeric(LTR$Trait)
LTR <- LTR %>% 
  mutate(
    Position_bin = cut(Position,breaks = 500)
  )
ggplot(LTR) + 
  geom_tile(aes(x = Position_bin, y = Chromosome, fill = Trait)) +
  scale_fill_gradient(low = "white", high = "red") + 
  theme_Publication()





# LTR methylation ---------------------------------------------------------
aver <- function(data){mean(data,na.rm = T,trim = 0.0025)}
Dark2 <- brewer.pal(8,'Dark2')
Paired <- brewer.pal(8,'Paired')
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
LTR.CG_methy <- plot_data('../data/Epi/CG.LTR.cent.matrix.gz','CG methylation')
LTR.CG_densi <- plot_data('../data/Epi/LTR.CG_density.matrix.gz','CG densities')

ggplot() + 
  geom_line(data=LTR.CG_methy,aes(x=pos,y = value),color = Paired[1]) + 
  geom_line(data = LTR.CG_densi, aes(x = pos, y = rescale(value,c(0.83,0.97))),color = Paired[2])+
  scale_y_continuous(breaks = pretty_breaks(n = 5), # 坐标刻度，五个刻度，六个区间
                     sec.axis = sec_axis(
                       ~rescale(.,c(300,700)), # 将坐标轴重新scale
                       name = " The densities of CG")) +
  scale_x_continuous(
    breaks = c(1,90,150,240),
    labels = c('-9K','TSS','TES','9K'),
    expand = c(0.005,0.005))  + theme_Publication() + 
  theme(
    axis.ticks.y.right = element_line(color = Paired[2],size = 1),
    axis.line.y.right = element_line(color = Paired[2],size = 1),
    axis.title.y.right = element_text(color = Paired[2],size = 16),
    axis.ticks.y.left = element_line(color = Paired[1],size = 1),
    axis.line.y.left = element_line(color = Paired[1],size = 1),
    axis.title.y.left = element_text(color = Paired[1],size = 16)
  )
ggsave('../figure/Epi/LTR.CG.metaplot.pdf',width = 4.3,height = 3.8)



LTR.CHG_methy <- plot_data('../data/Epi/CHG.LTR.cent.matrix.gz','CHG methylation')
LTR.CHG_densi <- plot_data('../data/Epi/LTR.CHG_density.matrix.gz','CHG densities')


ggplot() + 
  geom_line(data=LTR.CHG_methy,aes(x=pos,y = value),color = Paired[3]) + 
  geom_line(data = LTR.CHG_densi, aes(x = pos, y = rescale(value,c(0.47,0.7))),color = Paired[4])+
  scale_y_continuous(breaks = pretty_breaks(n = 5), # 坐标刻度，五个刻度，六个区间
                     sec.axis = sec_axis(
                       ~rescale(.,c(300,700)), # 将坐标轴重新scale
                       name = " The densities of CHG")) +
  scale_x_continuous(
    breaks = c(1,90,150,240),
    labels = c('-9K','TSS','TES','9K'),
    expand = c(0.005,0.005))  + theme_Publication() + 
  theme(
    axis.ticks.y.right = element_line(color = Paired[4],size = 1),
    axis.line.y.right = element_line(color = Paired[4],size = 1),
    axis.title.y.right = element_text(color = Paired[4],size = 16),
    axis.ticks.y.left = element_line(color = Paired[3],size = 1),
    axis.line.y.left = element_line(color = Paired[3],size = 1),
    axis.title.y.left = element_text(color = Paired[3],size = 16)
  )
ggsave('../figure/Epi/LTR.CHG.metaplot.pdf',width = 4.3,height = 3.8)





LTR.CHH_methy <- plot_data('../data/Epi/CHH.LTR.cent.matrix.gz','CHH methylation')
LTR.CHH_densi <- plot_data('../data/Epi/LTR.CHH_density.matrix.gz','CHH densities')


ggplot() + 
  geom_line(data=LTR.CHH_methy,aes(x=pos,y = value),color = Paired[5]) + 
  geom_line(data = LTR.CHH_densi, aes(x = pos, y = rescale(value,c(0.014,0.028))),color = Paired[6])+
  scale_y_continuous(breaks = pretty_breaks(n = 5), # 坐标刻度，五个刻度，六个区间
                     sec.axis = sec_axis(
                       ~rescale(.,c(300,700)), # 将坐标轴重新scale
                       name = " The densities of CHH")) +
  scale_x_continuous(
    breaks = c(1,90,150,240),
    labels = c('-9K','TSS','TES','9K'),
    expand = c(0.005,0.005))  + theme_Publication() + 
  theme(
    axis.ticks.y.right = element_line(color = Paired[6],size = 1),
    axis.line.y.right = element_line(color = Paired[6],size = 1),
    axis.title.y.right = element_text(color = Paired[6],size = 16),
    axis.ticks.y.left = element_line(color = Paired[5],size = 1),
    axis.line.y.left = element_line(color = Paired[5],size = 1),
    axis.title.y.left = element_text(color = Paired[5],size = 16)
  )
ggsave('../figure/Epi/LTR.CHH.metaplot.pdf',width = 4.3,height = 3.8)






# LTR methylation  vs. insertion time  ------------------------------------

# LTR.old.CG <- read_delim('../data/Epi/CG.LTR.intact.cent.low.bed',delim = '\t',col_names = F) %>% 
#   select(X1,X2,X3,X7) %>% mutate(Time = 'Old') %>% filter(X7!='.')
# LTR.mid.CG <- read_delim('../data/Epi/CG.LTR.intact.cent.mid.bed',delim = '\t',col_names = F) %>% 
#   select(X1,X2,X3,X7) %>% mutate(Time = 'Mid') %>% filter(X7!='.')
# LTR.young.CG <- read_delim('../data/Epi/CG.LTR.intact.cent.high.bed',delim = '\t',col_names = F) %>% 
#   select(X1,X2,X3,X7) %>% mutate(Time = 'Young') %>% filter(X7!='.')
# 
# LTR.CG <- rbind(LTR.old.CG,LTR.mid.CG,LTR.young.CG) %>% na.omit() %>% 
#   mutate(X7 = as.numeric(X7))
# 
# LTR.CG$Time <- factor(LTR.CG$Time,levels = c('Old','Mid','Young'))
# ggplot(LTR.CG,aes(x = Time, y = X7)) + geom_violin(aes(x = Time, y = X7,fill = Time, color = Time)) + 
#   geom_boxplot(aes(x = Time,y = X7,fill = Time),width = 0.3,outlier.alpha = 0) + 
#   theme_Publication() + 
#   geom_signif(comparisons=list(c("Mid", "Old")),map_signif_level = T,size = 0.9)+ 
#   geom_signif(comparisons=list(c("Mid", "Young")),map_signif_level = T,size = 0.9)
#   
# methyRatio.bar <- data.frame(
#   Type = c('Old','Mid','Young'),
#   CG = c(
#     mean(as.numeric(na.omit(LTR.old.CG$X7))),
#     mean(as.numeric(na.omit(LTR.mid.CG$X7))),
#     mean(as.numeric(na.omit(LTR.young.CG$X7)))
#   ),
#   CHH = c(
#     mean(as.numeric(na.omit(LTR.old.CHH$X7))),
#     mean(as.numeric(na.omit(LTR.mid.CHH$X7))),
#     mean(as.numeric(na.omit(LTR.young.CHH$X7)))
#   )
# )
# 
# 
# 
# LTR.old.CHH <- read_delim('../data/Epi/CHH.LTR.intact.cent.low.bed',delim = '\t',col_names = F) %>% 
#   select(X1,X2,X3,X7) %>% mutate(Time = 'Old') %>% filter(X7!='.')
# LTR.mid.CHH <- read_delim('../data/Epi/CHH.LTR.intact.cent.mid.bed',delim = '\t',col_names = F) %>% 
#   select(X1,X2,X3,X7) %>% mutate(Time = 'Mid') %>% filter(X7!='.')
# LTR.young.CHH <- read_delim('../data/Epi/CHH.LTR.intact.cent.high.bed',delim = '\t',col_names = F) %>% 
#   select(X1,X2,X3,X7) %>% mutate(Time = 'Young') %>% filter(X7!='.')
# 
# LTR.CHH <- rbind(LTR.old.CHH,LTR.old.CHH,LTR.young.CHH) %>% na.omit() %>% 
#   mutate(X7 = as.numeric(X7)) 
# LTR.CHH$Time <- factor(LTR.CHH$Time,levels = c('Old','Mid','Young'))
# ggplot(LTR.CHH,aes(x = Time, y = X7)) + geom_violin(aes(x = Time, y = X7,fill = Time, color = Time)) + 
#   geom_boxplot(aes(x = Time,y = X7,fill = Time),width = 0.3,outlier.alpha = 0) + 
#   theme_Publication() + 
#   geom_signif(comparisons=list(c("Mid", "Old")),map_signif_level = F,size = 0.9,y_position = 0.15)+ 
#   geom_signif(comparisons=list(c("Mid", "Young")),map_signif_level = T,size = 0.9)
# 
# mean(as.numeric(na.omit(LTR.old.CHH$X7)))/mean(as.numeric(na.omit(LTR.mid.CHH$X7)))





LTR.CG.same <- read_delim('../data/Epi/CG.LTR_1k.intact.cent.LTR_same.bed',delim = '\t',col_names = F) %>% 
  select(X1,X2,X3,X7) %>% mutate(Time = 'Same') %>% filter(X7!='.')
LTR.CG.diff <- read_delim('../data/Epi/CG.LTR_1k.intact.cent.LTR_diff.bed',delim = '\t',col_names = F) %>% 
  select(X1,X2,X3,X7) %>% mutate(Time = 'Diff') %>% filter(X7!='.')
LTR.CG.arm <- read_delim('../data/Epi/CG.LTR_1k.intact.arm.bed',delim = '\t',col_names = F) %>% 
  select(X1,X2,X3,X7) %>% mutate(Time = 'Arm') %>% filter(X7!='.')
LTR.CG <- rbind(LTR.CG.same,LTR.CG.diff,LTR.CG.arm) %>% na.omit() %>% 
  mutate(X7 = as.numeric(X7)) 

ggplot(LTR.CG) + stat_ecdf(aes(color = Time,x = X7)) + 
  xlim(c(0.8,1)) + theme_Publication() + scale_colour_Publication() + 
  labs(x = 'DNA methylation',y = 'Cumulative fraction')
summary(aov(X7~Time,data = LTR.CG))
ggsave('../figure/Epi/LTR.CG.ecdf.pdf',width = 3,height = 4)


LTR.CHH.same <- read_delim('../data/Epi/CHH.LTR_1k.intact.cent.LTR_same.bed',delim = '\t',col_names = F) %>% 
  select(X1,X2,X3,X7) %>% mutate(Time = 'Same') %>% filter(X7!='.')
LTR.CHH.diff <- read_delim('../data/Epi/CHH.LTR_1k.intact.cent.LTR_diff.bed',delim = '\t',col_names = F) %>% 
  select(X1,X2,X3,X7) %>% mutate(Time = 'Diff') %>% filter(X7!='.')
LTR.CHH.arm <- read_delim('../data/Epi/CHH.LTR_1k.intact.arm.bed',delim = '\t',col_names = F) %>% 
  select(X1,X2,X3,X7) %>% mutate(Time = 'Arm') %>% filter(X7!='.')
LTR.CHH <- rbind(LTR.CHH.same,LTR.CHH.diff,LTR.CHH.arm) %>% na.omit() %>% 
  mutate(X7 = as.numeric(X7)) 

ggplot(LTR.CHH) + stat_ecdf(aes(color = Time,x = X7)) + 
  xlim(c(0,0.05)) + theme_Publication() + scale_colour_Publication() + 
  labs(x = 'DNA methylation',y = 'Cumulative fraction')

summary(aov(X7~Time,data = LTR.CHH))
ggsave('../figure/Epi/LTR.CHH.ecdf.pdf',width = 3,height = 4)


