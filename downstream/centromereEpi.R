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
# centromere to proximal -----------------------------------------------------
K4.cent <- plot_data('../data/Epi/cent.mid.H3K4me3.matrix.gz','H3K4me3') 
K9.cent <- plot_data('../data/Epi/cent.mid.H3K9me3.matrix.gz','H3K9me3')
K27.cent <- plot_data('../data/Epi/cent.mid.H3K27me3.matrix.gz','H3K27me3')
silence.cent <- rbind(K9.cent,K27.cent)


ggplot() + 
  geom_line(data=silence.cent,aes(x=pos,y = value,color = variable)) +
  scale_color_manual(values = Paired[1:2]) + 
  geom_line(data=K4.cent,aes(x = pos,y = rescale(value,c(0,2.5))),color = Set2[2]) + 
  scale_y_continuous(breaks = pretty_breaks(n = 5), # 坐标刻度，五个刻度，六个区间
                     sec.axis = sec_axis(
                       ~rescale(.,c(0.14,0.34)), # 将坐标轴重新scale
                       name = "Second y-axis")) + theme_Publication() + 
  theme(
    axis.ticks.y.right = element_line(color = Set2[2],size = 1),
    axis.line.y.right = element_line(color = Set2[2],size = 1),
    axis.title.y.right = element_text(color = Set2[2],size = 16),
    axis.ticks.y.left = element_line(color = Set1[2],size = 1),
    axis.line.y.left = element_line(color = Set1[2],size = 1),
    axis.title.y.left = element_text(color = Set1[2],size = 16)
  ) +
  scale_x_continuous(
    breaks = c(1,200,400),
    labels = c('-200Mb','Centromere','200Mb'),
    expand = c(0.005,0.005))  + 
  labs(x = 'Position (bp)',y = 'Mean relative level') 

ggsave('../figure/CENH3/CENH3.proximal.histone.metaplot.pdf',width = 5,height = 4)



CG.cent <- plot_data('../data/Epi/cent.mid.CG.matrix.gz','CG') 
CHG.cent <- plot_data('../data/Epi/cent.mid.CHG.matrix.gz','CHG')
CHH.cent <- plot_data('../data/Epi/cent.mid.CHH.matrix.gz','CHH')


Spectral <- brewer.pal(9,'Spectral')
ggplot() + 
  geom_line(data=CG.cent,aes(x=pos,y = value),color = Spectral[2])  + 
  theme_Publication()  +
  scale_x_continuous(
    breaks = c(1,200,400),
    labels = c('-200Mb','Centromere','200Mb'),
    expand = c(0.005,0.005))  + 
  labs(x = 'Position (bp)',y = 'Mean relative level') 
ggsave('../figure/CENH3/CENH3.proximal.CG.metaplot.pdf',width = 4,height = 3)



ggplot() + 
  geom_line(data=CHG.cent,aes(x=pos,y = value),color = Spectral[3])  + 
  theme_Publication()  +
  scale_x_continuous(
    breaks = c(1,200,400),
    labels = c('-200Mb','Centromere','200Mb'),
    expand = c(0.005,0.005))  + 
  labs(x = 'Position (bp)',y = 'Mean relative level') 
ggsave('../figure/CENH3/CENH3.proximal.CHG.metaplot.pdf',width = 4,height = 3)

ggplot() + 
  geom_line(data=CHH.cent,aes(x=pos,y = value),color = Spectral[4])  + 
  theme_Publication()  +
  scale_x_continuous(
    breaks = c(1,200,400),
    labels = c('-200Mb','Centromere','200Mb'),
    expand = c(0.005,0.005))  + 
  labs(x = 'Position (bp)',y = 'Mean relative level') 
ggsave('../figure/CENH3/CENH3.proximal.CHH.metaplot.pdf',width = 4,height = 3)

CG.cent <- plot_data('../data/Epi/cent_CG.20M.matrix.gz','CG') 
CHG.cent <- plot_data('../data/Epi/cent_CHG.20M.matrix.gz','CHG')
CHH.cent <- plot_data('../data/Epi/cent_CHH.20M.matrix.gz','CHH')


Spectral <- brewer.pal(9,'Spectral')
ggplot() + 
  geom_line(data=CG.cent,aes(x=pos,y = value),color = Spectral[2])  + 
  theme_Publication()  +
  scale_x_continuous(
    breaks = c(1,200,400),
    labels = c('-20Mb','Centromere','20Mb'),
    expand = c(0.005,0.005))  + 
  labs(x = 'Position (bp)',y = 'Mean relative level') + ylim(0.92, 0.975)
ggsave('../figure/CENH3/CENH3.proximal.CG.metaplot.20M.pdf',width = 4,height = 3)



ggplot() + 
  geom_line(data=CHG.cent,aes(x=pos,y = value),color = Spectral[3])  + 
  theme_Publication()  +
  scale_x_continuous(
    breaks = c(1,200,400),
    labels = c('-20Mb','Centromere','20Mb'),
    expand = c(0.005,0.005))  + 
  labs(x = 'Position (bp)',y = 'Mean relative level') 
ggsave('../figure/CENH3/CENH3.proximal.CHG.metaplot.20M.pdf',width = 4,height = 3)

ggplot() + 
  geom_line(data=CHH.cent,aes(x=pos,y = value),color = Spectral[4])  + 
  theme_Publication()  +
  scale_x_continuous(
    breaks = c(1,200,400),
    labels = c('-20Mb','Centromere','20Mb'),
    expand = c(0.005,0.005))  + 
  labs(x = 'Position (bp)',y = 'Mean relative level') 
ggsave('../figure/CENH3/CENH3.proximal.CHH.metaplot.20M.pdf',width = 4,height = 3)


# centromere CENH3--------------------------------------------------------------

CHH.cenh3  <- plot_data('../data/CENH3/CHH.cent_nuclear.matrix.gz','CHH')
CG.cenh3  <- plot_data('../data/CENH3/CG.cent_nuclear.matrix.gz','CG')
CHG.cenh3  <- plot_data('../data/CENH3/CHG.cent_nuclear.matrix.gz','CHG')


ggplot() + 
  geom_line(data=CG.cenh3,aes(x=pos,y = value),color = Dark2[1]) +
  scale_x_continuous(
    breaks = c(1,200,400),
    labels = c('-2K','TSS','2K'),
    expand = c(0.005,0.005)) + theme_Publication() + 
  labs(x = 'Position (bp)',y = 'Mean relative level') 
ggsave('../figure/CENH3/CG.cent_nuclear.metaplot.pdf',width = 3,height = 3)

ggplot() + 
  geom_line(data=CHG.cenh3,aes(x=pos,y = value),color = Dark2[2]) +
  scale_x_continuous(
    breaks = c(1,100,200),
    labels = c('-2K','TSS','2K'),
    expand = c(0.005,0.005)) + theme_Publication()+ 
  labs(x = 'Position (bp)',y = 'Mean relative level') 
ggsave('../figure/CENH3/CHG.cent_nuclear.metaplot.pdf',width = 3,height = 3)

ggplot() + 
  geom_line(data=CHH.cenh3,aes(x=pos,y = value),color = Dark2[3]) +
  scale_x_continuous(
    breaks = c(1,100,200),
    labels = c('-2K','TSS','2K'),
    expand = c(0.005,0.005)) + theme_Publication()+ 
  labs(x = 'Position (bp)',y = 'Mean relative level') 
ggsave('../figure/CENH3/CHH.cent_nuclear.metaplot.pdf',width = 3,height = 3)

# 
# set.seed(123)
# CG.nuc <- read_delim('../data/CENH3/wn.CENH_nuclear.CG.bed',delim = '\t',col_names = F) %>% 
#   filter(X11!='.') %>% mutate(Ratio = as.numeric(X11),Type = 'Nucleosome') %>% select(Type,Ratio) 
# CG.nuc.sample <- CG.nuc[sample(1:dim(CG.nuc)[[1]],5000),]
# 
# 
# CG.nuc_flank <- read_delim('../data/CENH3/wn.CENH_nuclear.flank_2k.CG.bed',delim = '\t',col_names = F) %>% 
#   filter(X11!='.') %>% mutate(Ratio = X11,Type = 'Flank') %>% select(Type,Ratio)
# CG.nuc_flank.sample <- CG.nuc_flank[sample(1:dim(CG.nuc_flank)[[1]],5000),]
# CG <- rbind(CG.nuc.sample,CG.nuc_flank.sample) %>% mutate(Ratio = as.numeric(Ratio))
# 
# 
# 
# 
# CHG.nuc <- read_delim('../data/CENH3/wn.CENH_nuclear.CHG.bed',delim = '\t',col_names = F) %>% 
#   filter(X11!='.') %>% mutate(Ratio = as.numeric(X11),Type = 'Nucleosome') %>% select(Type,Ratio)
# CHG.nuc_flank <- read_delim('../data/CENH3/wn.CENH_nuclear.flank_2k.CHG.bed',delim = '\t',col_names = F) %>% 
#   filter(X11!='.') %>% mutate(Ratio = X11,Type = 'Flank') %>% select(Type,Ratio)
# CHG <- rbind(CHG.nuc,CHG.nuc_flank) %>% mutate(Ratio = as.numeric(Ratio))
# 
# CHH.nuc <- read_delim('../data/CENH3/wn.CENH_nuclear.CHH.bed',delim = '\t',col_names = F) %>% 
#   filter(X11!='.') %>% mutate(Ratio = as.numeric(X11),Type = 'Nucleosome') %>% select(Type,Ratio)
# CHH.nuc_flank <- read_delim('../data/CENH3/wn.CENH_nuclear.flank_2k.CHH.bed',delim = '\t',col_names = F) %>% 
#   filter(X11!='.') %>% mutate(Ratio = X11,Type = 'Flank') %>% select(Type,Ratio)
# CHH <- rbind(CHH.nuc,CHH.nuc_flank) %>% mutate(Ratio = as.numeric(Ratio))
# 
# 
# ggplot(CG) + stat_ecdf(aes(x = Ratio, color = Type)) + xlim(c(0.75,1))
# ggplot(CHG) + stat_ecdf(aes(x = Ratio, color = Type)) + xlim(c(0.75,1))
# ggplot(CHH) + stat_ecdf(aes(x = Ratio, color = Type)) + xlim(c(0,0.05))



# CENH3 region methylation ------------------------------------------------

CG.e <- read_delim('../data/CENH3/CENH3_depleted.CG.bg',delim = '\t',col_names = F) %>% 
  mutate(Type = 'Enriched',
         X4 = as.numeric(replace(X4, X4 == ".", 0)),
         X4 = replace_na(X4, 0))
CG.d <- read_delim('../data/CENH3/CENH3_enriched.CG.bg',delim = '\t',col_names = F) %>% 
  mutate(Type = 'Depleted',
         X4 = as.numeric(replace(X4, X4 == ".", 0)),
         X4 = replace_na(X4, 0))

CG <- rbind(CG.e,CG.d)
ggplot(CG) + stat_ecdf(aes(x = X4,color = Type))  + theme_Publication() + 
  xlim(c(0.9,0.98)) + 
  labs(x = 'Methylation level', y = 'Cumulative Fration')
ggsave('../figure/CENH3/CENH3.subdomain.ecdf.CG.pdf',width = 3,height = 3)

ggplot(CG,aes(x = Type,y = X4,fill = Type)) + geom_violin()  + ylim(c(0.8,1)) + 
  geom_boxplot(outlier.alpha = 0,width = 0.4) + theme_Publication() + 
  scale_fill_Publication('Accent') + 
  labs(x='Subdomain',y = 'Methylation level') + 
  theme(
    legend.position = 'None'
  ) + 
  geom_signif(comparisons=list(c("Depleted", "Enriched")),map_signif_level = T,size = 0.9,y_position = 0.8)
ggsave('../figure/CENH3/CENH3.subdomain.boxplot.CG.pdf',width = 2.5,height = 4)




CHG.e <- read_delim('../data/CENH3/CENH3_depleted.CHG.bg',delim = '\t',col_names = F) %>% 
  mutate(Type = 'Enriched',
         X4 = as.numeric(replace(X4, X4 == ".", 0)),
         X4 = replace_na(X4, 0))
CHG.d <- read_delim('../data/CENH3/CENH3_enriched.CHG.bg',delim = '\t',col_names = F) %>% 
  mutate(Type = 'Depleted',
         X4 = as.numeric(replace(X4, X4 == ".", 0)),
         X4 = replace_na(X4, 0))

CHG <- rbind(CHG.e,CHG.d)
ggplot(CHG) + stat_ecdf(aes(x = X4,color = Type))  + theme_Publication() + 
  xlim(c(0.3,0.8)) + 
  labs(x = 'Methylation level', y = 'Cumulative Fration')
ggsave('../figure/CENH3/CENH3.subdomain.ecdf.CHG.pdf',width = 3,height = 3)


ggplot(CHG,aes(x = Type,y = X4,fill = Type)) + geom_violin()  + 
  geom_boxplot(outlier.alpha = 0,width = 0.4) + theme_Publication() + 
  scale_fill_Publication('Accent') + 
  labs(x='Subdomain',y = 'Methylation level') + 
  theme(
    legend.position = 'None'
  ) + 
  geom_signif(comparisons=list(c("Depleted", "Enriched")),map_signif_level = T,size = 0.9,y_position = 0.92)
ggsave('../figure/CENH3/CENH3.subdomain.boxplot.CHG.pdf',width = 2.5,height = 4)



CHH.e <- read_delim('../data/CENH3/CENH3_depleted.CHH.bg',delim = '\t',col_names = F) %>% 
  mutate(Type = 'Enriched',
         X4 = as.numeric(replace(X4, X4 == ".", 0)),
         X4 = replace_na(X4, 0))
CHH.d <- read_delim('../data/CENH3/CENH3_enriched.CHH.bg',delim = '\t',col_names = F) %>% 
  mutate(Type = 'Depleted',
         X4 = as.numeric(replace(X4, X4 == ".", 0)),
         X4 = replace_na(X4, 0))

CHH <- rbind(CHH.e,CHH.d)
ggplot(CHH) + stat_ecdf(aes(x = X4,color = Type)) + xlim(c(0.005,0.03)) + theme_Publication() + 
  labs(x = 'Methylation level', y = 'Cumulative Fration')

ggsave('../figure/CENH3/CENH3.subdomain.ecdf.CHH.pdf',width = 3,height = 3)


ggplot(CHH,aes(x = Type,y = X4,fill = Type)) + geom_violin()  + ylim(c(0,0.05)) + 
  geom_boxplot(outlier.alpha = 0,width = 0.4) + theme_Publication() + 
  scale_fill_Publication('Accent') + 
  labs(x='Subdomain',y = 'Methylation level') + 
  theme(
    legend.position = 'None'
  ) + 
  geom_signif(comparisons=list(c("Depleted", "Enriched")),map_signif_level = T,size = 0.9,y_position = 0.03)
ggsave('../figure/CENH3/CENH3.subdomain.boxplot.CHH.pdf',width = 2.5,height = 4)


# subdomains --------------------------------------------------------------
CENH3 <- read_delim('../data/CENH3/wn_CENH3.100k.cent.bg',delim = '\t',col_names = F)
Input <- read_delim('../data/CENH3/wn_input.100k.cent.bg',delim = '\t',col_names = F)
mean(Input$X4) # 0.04779983
sd(Input$X4) # 0.02635133

CENH3$p <- apply(CENH3[,4], 1, function(x){pnorm(x,mean = 0.04779983,sd = 0.02635133,lower.tail = F)})
write_delim(CENH3,'../data/CENH3/CENH3_enrich_pvalue.100k.bed',delim = '\t',col_names = F)




CG <- read_delim('../data/CENH3/CENH3_enrich_pvalue.CG.bed',delim = '\t',col_names = F) %>% 
  filter(X6!='.') %>% mutate(X4 = as.numeric(X5),X6 = as.numeric(X6))
ggplot(CG) + geom_point(aes(x =log10(X5),y=(X6)) ) + ylim(c(0.75,1))


lm(X5~X6,data = CG)




# CRM epi -----------------------------------------------------------------

CRM.cenh3 <- plot_data('../data/CENH3/CRM.cenh3.matrix.gz','CRM') 
random.cenh3 <- plot_data('../data/CENH3/CRM_random_seed2.cenh3.matrix.gz','Random') 


CENH3.signal <- rbind(CRM.cenh3,random.cenh3)
ggplot(CENH3.signal) + geom_line(aes(x = pos,y = value, color = variable),size = 0.8) + 
  theme_Publication() + scale_colour_Publication('Set1') + 
  labs(x = 'Position (bp)',
       y = 'Mean relative signal\n',
       color = '')  +
  scale_x_continuous(
    breaks = c(1,60,90,150),
    labels = c('-6K','TSS','TES','6K'),
    expand = c(0.005,0.005)) + 
  theme(
    legend.position = c(0.8,1),
    legend.direction = 'vertical'
  )
ggsave('../figure/CENH3/CRM.cenh3.metaplot.pdf',width = 4,height = 3)





CRM.CG <- read_delim('../data/CENH3/CRM.cent.CG.bg',delim = '\t',col_names = F) %>% 
  mutate(Type = 'CRM') 
random.CG <- read_delim('../data/CENH3/CRM_random_seed2.cent.CG.bg',delim = '\t',col_names = F) %>% 
  mutate(Type = 'Random') 

CG.signal <- rbind(CRM.CG,random.CG) %>% filter(X7!=".") %>% mutate(X7 = as.numeric(X7))
ggplot(CG.signal,aes(x = X7, fill = Type,color = Type)) + stat_ecdf() + xlim(0.75,0.99) + 
  theme_Publication() + scale_colour_Publication('Accent')




ggplot(CG.signal,aes(x = X7, fill = Type,color = Type)) + geom_density(alpha = 0.3) + 
  scale_fill_Publication('Set2') + theme_Publication() + 
  scale_colour_Publication('Set2') + labs(x = 'Methylation level', y = 'Relative density')
ggsave('../figure/CENH3/CRM.CG.density.pdf',width = 4.3,height = 3.5)
ggplot(CG.signal,aes(x = X7, fill = Type,color = Type),) + geom_density(alpha = 0.3) + 
  scale_fill_Publication('Set2') + theme_Publication() + 
  coord_cartesian(xlim=c(0.8,1.0)) + scale_colour_Publication('Set2') + 
  labs(x = 'Methylation level', y = 'Relative density')
ggsave('../figure/CENH3/CRM.CG.local.density.pdf',width = 4.3,height = 3.5)


CRM.CHG <- read_delim('../data/CENH3/CRM.cent.CHG.bg',delim = '\t',col_names = F) %>% 
  mutate(Type = 'CRM') 
random.CHG <- read_delim('../data/CENH3/CRM_random_seed2.cent.CHG.bg',delim = '\t',col_names = F) %>% 
  mutate(Type = 'Random') 

CHG.signal <- rbind(CRM.CHG,random.CHG) %>% filter(X7!=".") %>% mutate(X7 = as.numeric(X7))
ggplot(CHG.signal,aes(x = X7, fill = Type,color = Type)) + geom_density(alpha = 0.3) + 
  scale_fill_Publication('Set2') + theme_Publication() + 
  labs(x = 'Methylation level', y = 'Density') + scale_colour_Publication()
ggsave('../figure/CENH3/CRM.CHG.density.pdf',width = 4.3,height = 3.5)

CRM.CHH <- read_delim('../data/CENH3/CRM.cent.CHH.bg',delim = '\t',col_names = F) %>% 
  mutate(Type = 'CRM') 
random.CHH <- read_delim('../data/CENH3/CRM_random_seed2.cent.CHH.bg',delim = '\t',col_names = F) %>% 
  mutate(Type = 'Random') 

CHH.signal <- rbind(CRM.CHH,random.CHH) %>% filter(X7!=".") %>% mutate(X7 = as.numeric(X7))
ggplot(CHH.signal,aes(x = X7, fill = Type,color = Type),) + geom_density(alpha = 0.3) + 
  scale_fill_Publication('Set2') + theme_Publication()+ 
  scale_colour_Publication('Set2') + labs(x = 'Methylation level', y = 'Relative density')
ggsave('../figure/CENH3/CRM.CHH.density.pdf',width = 4.3,height = 3.5)

ggplot(CHH.signal,aes(x = X7, fill = Type,color = Type),) + geom_density(alpha = 0.3) + 
  scale_fill_Publication('Set2') + theme_Publication() + 
  coord_cartesian(xlim=c(0.0,0.05))+ 
  scale_colour_Publication('Set2') + labs(x = 'Methylation level', y = 'Relative density')
ggsave('../figure/CENH3/CRM.CHH.local.density.pdf',width = 4.3,height = 3.5)



CRM.drip <- read_delim('../data/CENH3/CRM.cent.DRIP.bg',delim = '\t',col_names = F) %>% 
  mutate(Type = 'CRM') 
random.drip <- read_delim('../data/CENH3/CRM_random_seed2.cent.DRIP.bg',delim = '\t',col_names = F) %>% 
  mutate(Type = 'Random') 
drip.signal <- rbind(CRM.drip,random.drip) %>% filter(X7!=".") %>% mutate(X7 = as.numeric(X7))
ggplot(drip.signal,aes(x = Type, y = X7),) + geom_violin(aes(fill = Type)) + 
  geom_boxplot(width = 0.4,outlier.alpha = 0) + ylim(c(0,3)) + 
  theme_Publication() + theme(
    
  )


ggplot(drip.signal,aes(color = Type, x = X7),) + stat_ecdf() + xlim(0.1,13) + 
  theme_Publication() + scale_colour_Publication('Dark2') + 
  labs(x = 'Mean signal level', y = 'Cumulative fraction') + 
  theme(
    legend.position = c(0.7,0.3),
    legend.direction = 'vertical'
  )
ggsave('../figure/CENH3/CRM.Rloop.ecdf.pdf',width = 3.5,height = 3.5)


library(ggpubr) 
ggplot(drip.signal,aes(x = Type, y = X7),) + geom_violin(aes(fill = Type)) +  
  geom_boxplot(width = 0.15,outlier.alpha = 0)+  ylim(0.1,14) + theme_Publication() + 
  scale_colour_Publication() + 
  geom_signif(comparisons = list(c("CRM", "Random")), y_position = 12,
              map_signif_level=F) + 
  labs(y = 'Mean signal level')+ scale_fill_Publication('Dark2') 
ggsave('../figure/CENH3/CRM.Rloop.violin.pdf',width = 3,height = 4)
  

