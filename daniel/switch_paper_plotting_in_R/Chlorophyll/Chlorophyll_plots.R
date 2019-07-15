##don't forget the library dependencies
require(dplyr)
require(ggplot2)
library(ggpubr)
theme_set(theme_classic())

msr <- read.csv("/home/daniel/Sync/Google Drive/Niyogi Lab/Zofingiensis/Switch Paper/Switch Paper Part B/Replot/Chlorophyll_per_Cell.csv", header=TRUE,as.is=TRUE)
msr$strain = factor(msr$strain, levels= c("WT", "hxk1-1", "hxk1-2"))
msr1 <- msr %>% group_by(strain,glc) %>%  dplyr::summarise(average = mean(chl_per_cell), sd = sd(chl_per_cell))
msr$interact <-interaction(msr$strain,msr$glc)


p1<-  ggplot() +
  geom_bar(data=msr1, aes(x= interaction(glc,strain),y=average, fill=glc), color="black",position = "dodge", stat="identity", alpha = 0.8, width = 0.6)+scale_fill_manual(values=c("grey95", "grey15"))+
  geom_jitter(data=msr,aes(x= interaction(glc,strain),y=chl_per_cell),size = 2, color ="black", fill="white",stroke=0.5, alpha = 0.8, shape=21,width = 0.2)+
  geom_errorbar(data = msr1, aes(x=interaction(glc,strain),ymin=average- sd, ymax=average +sd),width = 0.2)
xlabels <-c("WT(-)", "WT(+)", expression(paste(italic("hxk1"),"-1(-)")), expression(paste(italic("hxk1"),"-1(+)")), expression(paste(italic("hxk1"),"-2(-)")),expression(paste(italic("hxk1"),"-2 (+)")))
p1 +scale_x_discrete(labels = xlabels)+ggtitle("Chlorophyll Per Cell")+ ylab(label = "fmol/cell") + xlab(label = " ")+ theme_classic(base_size = 20)+
  theme(plot.title = element_text(color = "black", size = 30, hjust = 0.5)) +scale_y_continuous(breaks = seq(0,0.7,0.10),limits = c(0,0.70),expand=c(0,0),)

#clustered 
p11<- ggplot(msr1, aes(x=strain,y=average,fill=glc,ymin=average-sd,ymax=average+sd))+
  geom_bar(position = "dodge",stat = "identity",color="black",alpha=0.8,width=0.5)+scale_fill_manual(values=c("grey95", "grey15"))+
  geom_point(data=msr,aes(x= strain,y=chl_per_cell,group=msr$interact),inherit.aes = FALSE,shape=21,color="black",fill="white",
             position = position_dodge(width=0.5),alpha = 0.7,size=4,stroke=1.1)+
  geom_errorbar(width=0.3, position=position_dodge(0.5),size=0.6)+scale_y_continuous(breaks = seq(0,0.7,0.1),limits = c(0,0.7),expand=c(0,0))+
  xlab(label = " ")+
  theme(legend.position = "none")+
  ggtitle("Chlorophyll Per Cell")
p11




## Total Chlorophyll

msr2 <- read.csv("/home/daniel/Sync/Google Drive/Niyogi Lab/Zofingiensis/Switch Paper/Switch Paper Part B/Replot/Chlorophyll_Total.csv", header=TRUE,as.is=TRUE)
msr2$strain = factor(msr2$strain, levels= c("WT", "hxk1-1", "hxk1-2"))
msr3 <- msr2 %>% group_by(strain,glc) %>%  dplyr::summarise(average = mean(total_chl), sd = sd(total_chl))
msr2$interact <-interaction(msr$strain,msr$glc)                                                            
                                                            
 
p12<- ggplot(msr3, aes(x=strain,y=average,fill=glc,ymin=average-sd,ymax=average+sd))+
  geom_bar(position = "dodge",stat = "identity",color="black",alpha=0.8,width=0.5)+scale_fill_manual(values=c("grey95", "grey15"))+
  geom_point(data=msr2,aes(x= strain,y=total_chl,group=msr2$interact),inherit.aes = FALSE,shape=21,color="black",fill="white",position = position_dodge(width=0.5),alpha = 0.7,size=4,stroke=1.1)+
  geom_errorbar(width=0.3, position=position_dodge(0.5),size=0.6)+scale_y_continuous(breaks = seq(0,12,2),limits = c(0,12),expand=c(0,0))+
  xlab(label = " ")+
  theme(legend.position = "none")+
  ggtitle("Total Chlorophyll")
p12
