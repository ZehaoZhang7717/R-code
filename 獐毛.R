library(readxl)
library(ggpubr)
library(dplyr)
dat <- data.frame(read_xlsx('新建文件夹/数据文件.xlsx',3))
rownames(dat) <-dat[,1] 
dat <- dat[,-1]
group <- data.frame(read_xlsx('新建文件夹/数据文件.xlsx',2))
rownames(group) <- group[,1] 
group <- group[,-1]
dat <- dat[rownames(group),]
data <- cbind(dat,group)
data$Treatment <- factor(data$Treatment,levels = c('N0P0','N1P1','N2P1','N3P1'))
####阿尔法多样性单因素方差分析多重比较
plant1 <- unique(data$Plant)%>% as.vector()
soil1 <- unique(data$Soil)%>% as.vector()
library(reshape2)
library(agricolae)
library(multcomp)
library(tidyverse)
result <- data.frame()
a <- 1

df <- data[data$Plant==plant1[a],]
for (b in 1:2) {
  df2 <- df[df$Soil==soil1[b],]
  df2 <- melt(df2,id=c('Plant','Treatment','Soil'))
  list <- unique(df2$variable) %>% as.vector() 
  for(i in 1:length(list)){
    df1<- df2 %>%
      dplyr::filter(variable == list[i])
    aov<- aov( value~Treatment, df1)
    lsd<- LSD.test(aov, "Treatment")
    lsdr<- lsd$groups %>% as.data.frame() %>%
      mutate(Treatment= rownames(lsd$groups),
             variable= list[i]) %>%
      dplyr::select(variable, Treatment, groups,value)
    lsdr$Soil <- soil1[b]
    lsdr$Plant <- plant1[a]
    result <- rbind(result,lsdr)
  }
}

library(ggplot2)
library(RColorBrewer)
library(lemon)
result$Treatment <- factor(result$Treatment,levels = c('N0P0','N1P1','N2P1','N3P1'))
data1 <- melt(data,id=c('Plant','Soil','Treatment'))
data1$Treatment  <- factor(data1$Treatment,levels = c('N0P0', 'N1P1','N2P1','N3P1'))

Chao1 <- ggplot(data1[data1$Plant==plant1[a]&data1$variable=='Chao1',],
                aes(Soil,value,color=Treatment,fill=Treatment))+
  stat_summary(fun = "mean", geom = "bar", width=0.5,
               position = position_dodge(0.7))+
  stat_summary(fun.data = "mean_sdl",fun.args = list(mult=1),
               geom = "errorbar",
               position = position_dodge(0.7),
               width = .2)+
  labs(y='Chao1')+
  theme_classic2()+
  theme(legend.title = element_blank(),
        axis.text.y = element_text(colour="black", size=13),
        axis.text.x = element_text(colour="black", size=13),
        axis.title = element_text(colour="black", size=16),
        legend.position = "bottom",
        legend.text = element_text(size=13),
        legend.key = element_rect(fill = NA),
        legend.key.size = unit(0.7, "cm"),
        legend.box.spacing = unit(0, "cm"),
        panel.background = element_blank(),
        plot.title = element_text(hjust=0.5, size=14),
        text=element_text(family="serif"))+
  geom_text(data = result[result$Plant==plant1[a]&result$variable=='Chao1',],aes(Soil,value,label=groups,color=Treatment),
            position = position_dodge(0.7),family='serif',vjust=0,
            size=6)

Shannon <- ggplot(data1[data1$Plant==plant1[a]&data1$variable=='Shannon',],
                  aes(Soil,value,color=Treatment,fill=Treatment))+
  stat_summary(fun = "mean", geom = "bar", width=0.5,
               position = position_dodge(0.7))+
  stat_summary(fun.data = "mean_sdl",fun.args = list(mult=1),
               geom = "errorbar",
               position = position_dodge(0.7),
               width = .2)+
  labs(y='Shannon')+
  theme_classic2()+
  theme(legend.title = element_blank(),
        axis.text.y = element_text(colour="black", size=13),
        axis.text.x = element_text(colour="black", size=13),
        axis.title = element_text(colour="black", size=16),
        legend.position = "bottom",
        legend.text = element_text(size=13),
        legend.key = element_rect(fill = NA),
        legend.key.size = unit(0.7, "cm"),
        legend.box.spacing = unit(0, "cm"),
        panel.background = element_blank(),
        plot.title = element_text(hjust=0.5, size=14),
        text=element_text(family="serif"))+
  geom_text(data = result[result$Plant==plant1[a]&result$variable=='Shannon',],aes(x = Soil,y = value,label=groups,color=Treatment),
            position = position_dodge(0.7),family='serif',vjust=0,
            size=6)
grid_arrange_shared_legend( Chao1, Shannon, 
                                   ncol = 2, nrow = 1,position='bottom')

##pcoa
library(readxl)
library(vegan)
otu <- data.frame(read_xlsx('新建文件夹/数据文件.xlsx',1))
rownames(otu) <- otu[,1]
otu <- otu[,-1]
group <- data.frame(read_xlsx('新建文件夹/数据文件.xlsx',2))
rownames(group) <- group[,1]
group <- group[,-1]
plant1 <- unique(group$Plant)
soil1 <- unique(group$Soil)
i=1

###根际与非根际分开 RS
dataa <- otu[,rownames(group[group$Plant==plant1[i]&group$Soil==soil1[1],])]
data1 <- vegdist(t(dataa),method = "bray")
pcoa <- cmdscale(data1, k=3, eig=T)
pcoa_points <- as.data.frame(pcoa$points)
sum_eig <- sum(pcoa$eig)
eig_percent <- round(pcoa$eig/sum_eig*100,1)
colnames(pcoa_points) <- paste0("PCoA", 1:3)
xlab=paste("PCoA1(",eig_percent[1],"%)", sep="")
ylab=paste("PCoA2(",eig_percent[2],"%)", sep="")
print(rownames(pcoa_points)==rownames(group[group$Plant==plant1[i]&group$Soil==soil1[1],]))
pcoa_result <- cbind(pcoa_points,group[group$Plant==plant1[i]&group$Soil==soil1[1],])
pcoa_result$Treatment <- factor(pcoa_result$Treatment,levels = c('N0P0',
                                                                 'N0P1','N1P1','N2P1','N3P1'))
set.seed(100)
dune.div <-adonis2(t(dataa) ~ Treatment, data = group[group$Soil=='RS',], permutations = 999, method="bray")
a <- group[group$Soil=='RS',]
dune.div2 <- anosim(t(dataa),a$Treatment ,permutations = 999,distance = 'bray')

library(ggpubr)
library(plyr)
library(gglayer)

str(pcoa_result)

RS <- ggplot(pcoa_result,aes(PCoA1,PCoA2))+
  geom_hline(yintercept=0,linetype = 3,size = 1) +
  geom_vline(xintercept=0,linetype = 3,size = 1)+
  geom_point(aes(fill=Treatment),size=3,shape=22)+
  geom_polygon(aes(fill = Treatment),color = 'black', alpha = 0.75,size=0.5, show.legend = FALSE)+
  labs(x=xlab,y=ylab,title = 'RS')+
  theme_test()+
  theme(legend.title = element_blank(),
        axis.text.y = element_text(colour="black", size=13),
        axis.text.x = element_text(colour="black", size=13),
        axis.title = element_text(colour="black", size=16),
        legend.position = "bottom",
        legend.text = element_text(size=13),
        legend.key = element_rect(fill = NA),
        legend.key.size = unit(0.7, "cm"),
        legend.box.spacing = unit(0, "cm"),
        panel.background = element_blank(),
        plot.title = element_text(hjust=0.5, size=14),
        text=element_text(family="serif"))
RS
##BS
dataa <- otu[,rownames(group[group$Plant==plant1[i]&group$Soil==soil1[2],])]
data1 <- vegdist(t(dataa),method = "bray")
pcoa <- cmdscale(data1, k=3, eig=T)
pcoa_points <- as.data.frame(pcoa$points)
sum_eig <- sum(pcoa$eig)
eig_percent <- round(pcoa$eig/sum_eig*100,1)
colnames(pcoa_points) <- paste0("PCoA", 1:3)
xlab=paste("PCoA1(",eig_percent[1],"%)", sep="")
ylab=paste("PCoA2(",eig_percent[2],"%)", sep="")
print(rownames(pcoa_points)==rownames(group[group$Plant==plant1[i]&group$Soil==soil1[2],]))
pcoa_result <- cbind(pcoa_points,group[group$Plant==plant1[i]&group$Soil==soil1[2],])
pcoa_result$Treatment <- factor(pcoa_result$Treatment,levels = c('N0P0',
                                                                 'N0P1','N1P1','N2P1','N3P1'))

set.seed(100)
dune.div <-adonis2(t(dataa) ~ Treatment, data = group[group$Soil=='BS',], permutations = 999, method="bray")
a <- group[group$Soil=='BS',]
dune.div2 <- anosim(t(dataa),a$Treatment ,permutations = 999,distance = 'bray')

library(ggpubr)
BS <- ggplot(pcoa_result,aes(PCoA1,PCoA2))+
  geom_hline(yintercept=0,linetype = 3,size = 1) +
  geom_vline(xintercept=0,linetype = 3,size = 1)+
  geom_point(aes(fill=Treatment),size=3,shape=22)+
  geom_polygon(aes(fill = Treatment),color = 'black', alpha = 0.75,size=0.5, show.legend = FALSE)+
  labs(x=xlab,y=ylab,title = 'BS')+
  theme_test()+
  theme(legend.title = element_blank(),
        axis.text.y = element_text(colour="black", size=13),
        axis.text.x = element_text(colour="black", size=13),
        axis.title = element_text(colour="black", size=16),
        legend.position = "bottom",
        legend.text = element_text(size=13),
        legend.key = element_rect(fill = NA),
        legend.key.size = unit(0.7, "cm"),
        legend.box.spacing = unit(0, "cm"),
        panel.background = element_blank(),
        plot.title = element_text(hjust=0.5, size=14),
        text=element_text(family="serif"))
library(lemon)
abc <- grid_arrange_shared_legend(RS,BS, ncol = 2, nrow = 1,position='bottom')
#ggsave(abc,filename = paste0(plant1[i],'.pdf'),width = 13,height = 6.5)






##pcoa
library(readxl)
library(vegan)
otu <- data.frame(read_xlsx('新建文件夹/数据文件.xlsx',1))
rownames(otu) <- otu[,1]
otu <- otu[,-1]
group <- data.frame(read_xlsx('新建文件夹/数据文件.xlsx',2))
rownames(group) <- group[,1]
group <- group[,-1]
plant1 <- unique(group$Plant)
soil1 <- unique(group$Soil)
i=1
###根际与非根际分开 RS
dataa <- otu[,rownames(group[group$Plant==plant1[i]&group$Soil=='RS',])]
data1 <- vegdist(t(dataa),method = "bray")
pcoa <- cmdscale(data1, k=3, eig=T)
pcoa_points <- as.data.frame(pcoa$points)
sum_eig <- sum(pcoa$eig)
eig_percent <- round(pcoa$eig/sum_eig*100,1)
colnames(pcoa_points) <- paste0("PCoA", 1:3)
xlab=paste("PCoA1(",eig_percent[1],"%)", sep="")
ylab=paste("PCoA2(",eig_percent[2],"%)", sep="")
print(rownames(pcoa_points)==rownames(group[group$Plant==plant1[i]&group$Soil==soil1[1],]))
pcoa_result <- cbind(pcoa_points,group[group$Plant==plant1[i]&group$Soil==soil1[1],])
pcoa_result$Treatment <- factor(pcoa_result$Treatment,levels = c('N0P0',
                                                                 'N1P1','N2P1','N3P1'))

plotinfo2 <- NULL
for (i in unique(pcoa_result$Treatment)) {
  tmp <- pcoa_result[pcoa_result$Treatment == i,] # 取出当前亚型当前来源下的数据
  avgx <- mean(tmp$PCoA1) # 计算横坐标均值
  avgy <- mean(tmp$PCoA2) # 计算纵坐标均值
  sdx <- sd(tmp$PCoA1) # 计算横坐标标准差
  sdy <- sd(tmp$PCoA2) # 计算纵坐标标准差
  
  plotinfo2 <- rbind.data.frame(plotinfo2,
                                data.frame(Treatment = i,
                                          
                                           avgx = avgx, # 添加圆的x位置
                                           avgy = avgy, # 添加圆的y位置
                                           sdx = sdx, # 添加圆的水平标准差
                                           sdy = sdy, # 添加圆的垂直标准差
                                           stringsAsFactors = F),
                                stringsAsFactors = F)}

library(ggplot2)
  RS <- ggplot(data = plotinfo2,aes(avgx,avgy))+
  geom_point(aes(fill=Treatment),shape=22,size=5)+
  geom_errorbar(aes(ymin = avgy-sdy,ymax = avgy+sdy,color=Treatment),width=0,show.legend = F) +
  geom_errorbarh(aes(xmin = avgx-sdx,xmax = avgx+sdx,color=Treatment),height=0,show.legend = F)+
  labs(x=xlab,y=ylab,title = 'RS')+
  theme_test()+
  theme(legend.title = element_blank(),
        axis.text.y = element_text(colour="black", size=13),
        axis.text.x = element_text(colour="black", size=13),
        axis.title = element_text(colour="black", size=16),
        legend.position = "bottom",
        legend.text = element_text(size=13),
        legend.key = element_rect(fill = NA),
        legend.key.size = unit(0.7, "cm"),
        legend.box.spacing = unit(0, "cm"),
        panel.background = element_blank(),
        plot.title = element_text(hjust=0.5, size=14),
        text=element_text(family="serif"))
##BS
  otu <- data.frame(read_xlsx('新建文件夹/数据文件.xlsx',1))
  rownames(otu) <- otu[,1]
  otu <- otu[,-1]
  group <- data.frame(read_xlsx('新建文件夹/数据文件.xlsx',2))
  rownames(group) <- group[,1]
  group <- group[,-1]
  plant1 <- unique(group$Plant)
  soil1 <- unique(group$Soil)
  i=1
  dataa <- otu[,rownames(group[group$Plant==plant1[i]&group$Soil=='BS',])]
  data1 <- vegdist(t(dataa),method = "bray")
  pcoa <- cmdscale(data1, k=3, eig=T)
  pcoa_points <- as.data.frame(pcoa$points)
  sum_eig <- sum(pcoa$eig)
  eig_percent <- round(pcoa$eig/sum_eig*100,1)
  colnames(pcoa_points) <- paste0("PCoA", 1:3)
  xlab=paste("PCoA1(",eig_percent[1],"%)", sep="")
  ylab=paste("PCoA2(",eig_percent[2],"%)", sep="")
  print(rownames(pcoa_points)==rownames(group[group$Plant==plant1[i]&group$Soil==soil1[1],]))
  pcoa_result <- cbind(pcoa_points,group[group$Plant==plant1[i]&group$Soil==soil1[1],])
  pcoa_result$Treatment <- factor(pcoa_result$Treatment,levels = c('N0P0',
                                                                   'N1P1','N2P1','N3P1'))
  
  plotinfo2 <- NULL
  for (i in unique(pcoa_result$Treatment)) {
    tmp <- pcoa_result[pcoa_result$Treatment == i,] # 取出当前亚型当前来源下的数据
    avgx <- mean(tmp$PCoA1) # 计算横坐标均值
    avgy <- mean(tmp$PCoA2) # 计算纵坐标均值
    sdx <- sd(tmp$PCoA1) # 计算横坐标标准差
    sdy <- sd(tmp$PCoA2) # 计算纵坐标标准差
    
    plotinfo2 <- rbind.data.frame(plotinfo2,
                                  data.frame(Treatment = i,
                                             
                                             avgx = avgx, # 添加圆的x位置
                                             avgy = avgy, # 添加圆的y位置
                                             sdx = sdx, # 添加圆的水平标准差
                                             sdy = sdy, # 添加圆的垂直标准差
                                             stringsAsFactors = F),
                                  stringsAsFactors = F)}
  
  library(ggplot2)
    BS <- ggplot(data = plotinfo2,aes(avgx,avgy))+
    geom_point(aes(fill=Treatment),shape=22,size=5)+
    geom_errorbar(aes(ymin = avgy-sdy,ymax = avgy+sdy,color=Treatment),width=0,show.legend = F) +
    geom_errorbarh(aes(xmin = avgx-sdx,xmax = avgx+sdx,color=Treatment),height=0,show.legend = F)+
    labs(x=xlab,y=ylab,title = 'BS')+
    theme_test()+
    theme(legend.title = element_blank(),
          axis.text.y = element_text(colour="black", size=13),
          axis.text.x = element_text(colour="black", size=13),
          axis.title = element_text(colour="black", size=16),
          legend.position = "bottom",
          legend.text = element_text(size=13),
          legend.key = element_rect(fill = NA),
          legend.key.size = unit(0.7, "cm"),
          legend.box.spacing = unit(0, "cm"),
          panel.background = element_blank(),
          plot.title = element_text(hjust=0.5, size=14),
          text=element_text(family="serif"))
    
BS
RS
library(lemon)
abc <- grid_arrange_shared_legend(RS,BS, ncol = 2, nrow = 1,position='bottom')
####门相对丰度表
library(readxl)
data <- data.frame(read_xlsx('新建文件夹/数据文件.xlsx',4))
rownames(data) <- data[,1]
data <- data[,-1]
data <- data.frame(t(data))
group <- data.frame(read_xlsx('新建文件夹/数据文件.xlsx',2))
rownames(group) <- group[,1]
group <- group[,-1]
library(reshape2)
data <- cbind(data,group)
dat <- melt(data,id=c('Plant','Soil','Treatment'))

                                                                                               
plant1 <- unique(dat$Plant)%>% as.vector()
soil1 <- unique(dat$Soil)%>% as.vector()
library(reshape2)
library(agricolae)
library(multcomp)
library(tidyverse)
result <- data.frame()
a <- 1

df <- dat[dat$Plant==plant1[1],]
for (b in 1:2) {
  df2 <- df[df$Soil==soil1[b],]
  list <- unique(df2$variable) %>% as.vector() 
  for(i in 1:9){
    df1<- df2 %>%
      dplyr::filter(variable == list[i])
    aov<- aov(value~Treatment, df1)
    lsd<- LSD.test(aov, "Treatment")
    lsdr<- lsd$groups %>% as.data.frame() %>%
      mutate(Treatment= rownames(lsd$groups),
             variable= list[i]) %>%
      dplyr::select(variable, Treatment, groups,value)
    lsdr$Soil <- soil1[b]
    lsdr$Plant <- plant1[a]
    result <- rbind(result,lsdr)
  }
}

sd <- aggregate(dat$value,by=list(dat$Soil,dat$Treatment,dat$variable),sd)
mean <- aggregate(dat$value,by=list(dat$Soil,dat$Treatment,dat$variable),mean)
write.csv(mean,file = 'mean.csv')
colnames(sd) <- c('Soil','Treatment','variable','sd')
data2 <- merge(result,sd,by = c('Soil','Treatment','variable'))

library(ggplot2)
str(dat)
dat$Treatment <- as.factor(dat$Treatment)
ggplot(dat[dat$variable!='Other',],
       aes(variable,value,color= Treatment,fill= Treatment))+
  stat_summary(fun = "mean", geom = "bar", width=0.5,
               position = position_dodge(0.7))+
  stat_summary(fun.data = 'mean_sdl',fun.args = list(mult=1),
               geom = "errorbar",
               position = position_dodge(0.7),
               width=0.2)+
  theme_test()+
  facet_rep_grid(Soil~.,scales = 'free_y')+
  geom_text(data = data2,aes(variable,y=value+sd,color=Treatment,label=groups),
            position = position_dodge(0.7),vjust=0)+
labs(x='Phylum',y = "Relative abundance (%)")+
theme(legend.title = element_blank(),
      axis.text.y = element_text(colour="black", size=13),
      axis.text.x = element_text(colour="black", size=13,angle = 45,
                                 vjust = 0.65),
      axis.title = element_text(colour="black", size=16),
      legend.position = "bottom",
      legend.text = element_text(size=13),
      legend.key = element_rect(fill = NA),
      legend.key.size = unit(0.7, "cm"),
      legend.box.spacing = unit(0, "cm"),
      panel.background = element_blank(),
      plot.title = element_text(hjust=0.4, size=14),
      text=element_text(family="serif"),
      strip.background = element_blank())
####距离衰减模型
library(readxl)
library(vegan)
otu <- data.frame(read_xlsx('新建文件夹/数据文件.xlsx',1))
rownames(otu) <- otu[,1]
otu <- otu[,-1]
group <- data.frame(read_xlsx('新建文件夹/数据文件.xlsx',2))
rownames(group) <- group[,1]
group <- group[,-1]
plant1 <- unique(group$Plant)
dis <- vegdist(t(otu[,rownames(group[group$Soil=='RS',])]), method = 'bray')
dis <- as.matrix(dis)
disR <- dis[,1:3]
disB <- vegdist(t(otu[,rownames(group[group$Soil=='BS',])]), method = 'bray')
disB <- as.matrix(disB)
disB <- disB[,1:3]

library(segmented)
library(reshape2)
disR <- melt(disR)
disB <- melt(disB)
DIS <- rbind(disR,disB)
colnames(DIS)[1] <- 'sample'
group$sample <- rownames(group)
dis1 <- merge(DIS,group,by = 'sample')
dis1$ratio <- rep(c(0,5,15,45),each=9,times=2)
dis1 <- dis1[dis1$value!=0,]
summary(lm(dis1[dis1$Soil=='BS',3]~dis1[dis1$Soil=='BS',7]))
dis2 <- dis1[dis1$Soil=='RS',]

fit_lm <- lm(value~ratio,dis2)
summary(fit_lm)
lm_seg1 <- segmented(fit_lm, seg.Z =~ratio, psi = 15)
library(ggplot2)
str(dis1)
ggplot(dis1,aes(ratio,value,fill=Soil))+
geom_smooth(aes(color=Soil),method = 'lm',formula = y~x,show.legend = F,alpha=0.2)+
  geom_point(size=5,shape=21,alpha=0.5)+
  scale_x_continuous(breaks=c(0,5,15,45))+
  theme_classic()+
  theme(legend.title = element_blank(),
        axis.text.y = element_text(colour="black", size=13),
        axis.text.x = element_text(colour="black", size=13),
        axis.title = element_text(colour="black", size=16),
        legend.position = "bottom",
        legend.text = element_text(size=13),
        legend.key = element_rect(fill = NA),
        legend.key.size = unit(0.7, "cm"),
        legend.box.spacing = unit(0, "cm"),
        panel.background = element_blank(),
        plot.title = element_text(hjust=0.5, size=14),
        text=element_text(family="serif"))+
  labs(y='Bray-Curtis dissimilarity',x='Nitrogen and phosphorus ratio')
  scale_color_manual(values = c("#f8756d",'#7cae00'))+
  scale_fill_manual(values = c("#f8756d",'#7cae00'))

###功能分析
library(readxl)
fun <- data.frame(read_xlsx('新建文件夹/数据文件.xlsx',9))
rownames(fun) <- fun[,1]
fun <- fun[,-1]
fun <- fun/colSums(fun)
fun <- data.frame(t(fun))
data <- fun

#write.csv(data,file = 'fun_relative_abundance.csv')


group <- data.frame(read_xlsx('新建文件夹/数据文件.xlsx',2))
rownames(group) <- group[,1]
group <- group[,-1]
library(reshape2)
data <- cbind(data,group)
dat <- melt(data,id=c('Plant','Soil','Treatment'))





plant1 <- unique(dat$Plant)%>% as.vector()
soil1 <- unique(dat$Soil)%>% as.vector()
library(reshape2)
library(agricolae)
library(multcomp)
library(tidyverse)
result <- data.frame()
a <- 1

df <- dat[dat$Plant==plant1[1],]
for (b in 1:2) {
  df2 <- df[df$Soil==soil1[b],]
  list <- unique(df2$variable) %>% as.vector() 
  for(i in 1:11){
    df1<- df2 %>%
      dplyr::filter(variable == list[i])
    aov<- aov(value~Treatment, df1)
    lsd<- LSD.test(aov, "Treatment")
    lsdr<- lsd$groups %>% as.data.frame() %>%
      mutate(Treatment= rownames(lsd$groups),
             variable= list[i]) %>%
      dplyr::select(variable, Treatment, groups,value)
    lsdr$Soil <- soil1[b]
    lsdr$Plant <- plant1[a]
    result <- rbind(result,lsdr)
  }
}

sd <- aggregate(dat$value,by=list(dat$Soil,dat$Treatment,dat$variable),sd)
mean <- aggregate(dat$value,by=list(dat$Soil,dat$Treatment,dat$variable),mean)
#write.csv(mean,file = 'mean_fun.csv')
colnames(sd) <- c('Soil','Treatment','variable','sd')
data2 <- merge(result,sd,by = c('Soil','Treatment','variable'))

library(ggplot2)
library(ggrepel)
library(lemon)
str(dat)
dat$Treatment <- as.factor(dat$Treatment)
ggplot(dat[dat$Soil=='BS',],
       aes(variable,value,fill= Treatment))+
  stat_summary(fun = "mean", geom = "bar", width=0.8,
               position = position_dodge(0.8),color= 'black')+
  stat_summary(fun.data = 'mean_sdl',fun.args = list(mult=1),
               geom = "errorbar",
               position = position_dodge(0.8),
               width=0.2)+
  ggplot2::coord_flip()+
  theme_test()+
  facet_rep_grid(Soil~.,scales = 'free_y')+
  geom_text(data = data2[data2$Soil=='BS',],aes(variable,y=value+sd,label=groups),size=4,
            position = position_dodge(0.8),hjust=-0.5,family='serif',vjust=0.2)+
  labs(x='Fun',y = "Relative abundance (%)")+
  theme(legend.title = element_blank(),
        axis.text.y = element_text(colour="black", size=13),
        axis.text.x = element_text(colour="black", size=13),
        axis.title = element_text(colour="black", size=16),
        legend.position = "bottom",
        legend.text = element_text(size=13),
        legend.key = element_rect(fill = NA),
        legend.key.size = unit(0.7, "cm"),
        legend.box.spacing = unit(0, "cm"),
        panel.background = element_blank(),
        plot.title = element_text(hjust=0.4, size=14),
        text=element_text(family="serif"),
        strip.background = element_blank())
aov<- aov(value~variable, dat)
lsd<- LSD.test(aov, "variable")


library(readxl)
library(tidyverse)
env <- data.frame(read_xlsx('新建文件夹/数据文件.xlsx',6))
group <- data.frame(read_xlsx('新建文件夹/数据文件.xlsx',2))
colnames(env)[1] <- 'sample'              
data <- merge(env,group,by = 'sample')
rownames(data) <- data[,1]
data <-  data[,-1]

data$Treatment <- factor(data$Treatment,levels = c('N0P0','N1P1','N2P1','N3P1'))

plant1 <- unique(data$Plant)%>% as.vector()
soil1 <- unique(data$Soil)%>% as.vector()
library(reshape2)
library(agricolae)
library(multcomp)
library(tidyverse)
result <- data.frame()
a <- 1

####属相对分析####
library(readxl)
data <- data.frame(read_xlsx('新建文件夹/数据文件.xlsx',5))
rownames(data) <- data[,1]
group <- data.frame(read_xlsx('新建文件夹/数据文件.xlsx',2))
data <- data.frame(t(data[,-1]))
dat <-cbind(group,data)
library(ggtree)
mean <-  aggregate(dat[,-1:-4],by = list(dat$Plant,dat$Treatment,dat$Soil),mean)
rownames(mean) <- paste0(mean$Group.3,mean$Group.2)
C <- mean[,-1:-3]%>%scale()%>%data.frame()
C <- C%>% mutate(B=row.names(.)) %>% melt()

 A <- mean[,-1:-3]%>% scale() %>% as.data.frame()
phr <- hclust(dist(A)) %>% 
  ggtree(layout="rectangular")
phc <- hclust(dist(t(A))) %>% 
  ggtree() 


library (ggplot2)
library (reshape2)#数据转换
require(scales)#数据缩放
library(ggtree)#聚类
library(aplot)#拼图
a <- phc[["data"]][["label"]]
b <-phc[["data"]][["y"]]
c <-phc[["data"]][["x"]]
datatt <- data.frame(a,b,c)
datatt <- datatt[order(datatt$b),]
AA <- na.omit(datatt)
  aaa <- AA$a
C$variable <- factor(C$variable,levels = aaa)
 ggplot(C,aes(B,variable,fill=value))+
  geom_raster()+
  scale_fill_gradient2(high="red", low="#c77cff", mid="white")+
  theme_bw()+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))+
  xlab(NULL) +
    ylab(NULL)+
  geom_vline(xintercept=4.5,size=0.5,linetype=1,color='black')+
   geom_hline(yintercept=7.5,size=0.5,linetype=1,color='black')+
   geom_hline(yintercept=15.5,size=0.5,linetype=1,color='black')+
   geom_hline(yintercept=18.5,size=0.5,linetype=1,color='black')+
   geom_hline(yintercept=28.5,size=0.5,linetype=1,color='black')+
   geom_hline(yintercept=28.5,size=0.5,linetype=1,color='black')+
   geom_hline(yintercept=40.5,size=0.5,linetype=1,color='black')+
   geom_hline(yintercept=37.5,size=0.5,linetype=1,color='black')+
   geom_hline(yintercept=45.5,size=0.5,linetype=1,color='black')+
   theme(legend.title = element_blank(),
         axis.text.y = element_text(colour="black", size=13),
         axis.text.x = element_text(colour="black", size=13),
         axis.title = element_text(colour="black", size=16),
         legend.text = element_text(size=13),
         legend.key = element_rect(fill = NA),
         legend.key.size = unit(0.7, "cm"),
         legend.box.spacing = unit(0, "cm"),
         panel.background = element_blank(),
         plot.title = element_text(hjust=0.4, size=14),
         text=element_text(family="serif"),
         strip.background = element_blank())

#属多因素方差分析
 library(readxl)
 data <- data.frame(read_xlsx('新建文件夹/数据文件.xlsx',5))
 rownames(data) <- data[,1]
 group <- data.frame(read_xlsx('新建文件夹/数据文件.xlsx',2))
 data <- data.frame(t(data[,-1]))
 dat <-cbind(group,data)
 dat$group <- paste0(dat$Soil,dat$Treatment)
library(reshape2) 
datlong <- melt(dat[,-1],id=c('Soil','Treatment','Plant','group')) 
soil1 <- unique(datlong$Soil) %>%as.vector()
Treatment <- unique(datlong$group) %>%as.vector()
list <- unique(datlong$variable)%>%as.vector()

result <- data.frame()
library(agricolae)
library(multcomp)
library(tidyverse)
for (s in 1:2) {
  dat_soil <- datlong[datlong$Soil==soil1[s],]
  for (i in 1:49) {
    dat_soil_var <- dat_soil[dat_soil$variable==list[i],]
    
    aov<- aov( value~group, dat_soil_var)
    lsd<- LSD.test(aov, "group")
    lsdr<- lsd$groups %>% as.data.frame() %>%
      mutate(group= rownames(lsd$groups),
             variable= list[i])
    lsdr$Soil <- soil1[s]
    result <- rbind(result,lsdr)
  }
}

result <- result[,c(2:4)]
colnames(result) <- c('group','B','variable')
D <- merge(C,result,by=c('B','variable'))
D$variable <- factor(D$variable,levels = aaa)
ggplot(D,aes(B,variable,fill=value))+
  geom_raster()+
  scale_fill_gradient2(high="red", low="#c77cff", mid="white")+
  theme_bw()+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))+
  xlab(NULL) +
  ylab(NULL)+
  geom_text(aes(B,variable,label=group),family='serif')+
  geom_vline(xintercept=4.5,size=0.5,linetype=1,color='black')+
  geom_hline(yintercept=7.5,size=0.5,linetype=1,color='black')+
  geom_hline(yintercept=15.5,size=0.5,linetype=1,color='black')+
  geom_hline(yintercept=18.5,size=0.5,linetype=1,color='black')+
  geom_hline(yintercept=28.5,size=0.5,linetype=1,color='black')+
  geom_hline(yintercept=28.5,size=0.5,linetype=1,color='black')+
  geom_hline(yintercept=40.5,size=0.5,linetype=1,color='black')+
  geom_hline(yintercept=37.5,size=0.5,linetype=1,color='black')+
  geom_hline(yintercept=45.5,size=0.5,linetype=1,color='black')+
  theme(legend.title = element_blank(),
        axis.text.y = element_text(colour="black", size=13),
        axis.text.x = element_text(colour="black", size=13),
        axis.title = element_text(colour="black", size=16),
        legend.text = element_text(size=13),
        legend.key = element_rect(fill = NA),
        legend.key.size = unit(0.7, "cm"),
        legend.box.spacing = unit(0, "cm"),
        panel.background = element_blank(),
        plot.title = element_text(hjust=0.4, size=14),
        text=element_text(family="serif"),
        strip.background = element_blank())

#相关性热图####
library(readxl)
df1 <- data.frame(read_xlsx('新建文件夹/数据文件.xlsx',6))
rownames(df1 ) <- df1[,1]
df1 <- df1[,-1]
group <- data.frame(read_xlsx('新建文件夹/数据文件.xlsx',2))
rownames(group) <- group[,1]
df2 <- data.frame(read_xlsx('新建文件夹/数据文件.xlsx',10))
rownames(df2) <- df2[,1]
df2 <- df2[,-1]

df1_R <- df1[rownames(group[group$Soil=='RS',]),]
df1_B <- df1[rownames(group[group$Soil=='BS',]),]

df2_R <- df2[rownames(group[group$Soil=='RS',]),]
df2_B <- df2[rownames(group[group$Soil=='BS',]),]

library(psych)
cor.result<-corr.test(df1_R,df2_R,method = "pearson")

library(tidyverse)
cor.result$p %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  pivot_longer(!rowname) %>% 
  mutate(p_value=case_when(
    value > 0.05 ~ "A",
    value >0.01 & value <= 0.05 ~ "B",
    value > 0.001 & value <= 0.01 ~ "D",
    value <= 0.001 ~ "E"
  )) -> new_df_R

cor.result$r %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  pivot_longer(!rowname) %>% 
  mutate(abs_cor=abs(value)) -> new_df2_R

cor.result<-corr.test(df1_B,df2_B,method = "pearson")

library(tidyverse)
cor.result$p %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  pivot_longer(!rowname) %>% 
  mutate(p_value=case_when(
    value > 0.05 ~ "A",
    value >0.01 & value <= 0.05 ~ "B",
    value > 0.001 & value <= 0.01 ~ "D",
    value <= 0.001 ~ "E"
  )) -> new_df_B
cor.result$r %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  pivot_longer(!rowname) %>% 
  mutate(abs_cor=abs(value)) -> new_df2_B

new_df_R$rowname <- paste0('R',new_df_R$rowname )
new_df_B$rowname <- paste0('B',new_df_B$rowname )


new_df2_R$rowname <- paste0('R',new_df2_R$rowname )
new_df2_B$rowname <- paste0('B',new_df2_B$rowname )

new_df1 <- rbind(new_df_R,new_df_B)
new_df2 <- rbind(new_df2_R,new_df2_B)
a <- colnames(df2)
new_df1$name <- factor(new_df1$name ,levels = a)
new_df2$name <- factor(new_df2$name ,levels = a)
new_df1$abs <-  abs(new_df1$value)
#write.csv(new_df2,file = 'corr.csv')
library(paletteer)
ggplot()+
  geom_tile(data=new_df1,
            aes(x=rowname,y=name,fill=p_value))+
  scale_fill_manual(values = c("white","#c0c0c0",
                               "#808080","#3f3f3f"))+

  scale_color_paletteer_c(palette = "ggthemes::Classic Red-Blue")+
  scale_alpha_manual(values = c(0,1,1,1))+
  scale_fill_manual(values = c("white","#c0c0c0",
                               "#808080","#3f3f3f"),
                    label=c(">0.05",
                            "0.01~0.05",
                            "0.001~0.01",
                            "<0.01"))+
  theme_bw()+
  theme(legend.key = element_rect(colour="black"),
        axis.text.x = element_text(angle = 45,hjust=0.5,vjust=0.5,
                                   colour="black", size=13),
        text = element_text(family = 'serif'),
        axis.text.y = element_text(colour="black", size=13),
        axis.title = element_text(colour="black", size=16),
        legend.text = element_text(size=13))+
geom_point(data=new_df2,
           aes(x=rowname,y=name,
               size=abs_cor,
               color=value))+
  geom_vline(xintercept=9.5,size=0.5,linetype=1,color='black')

