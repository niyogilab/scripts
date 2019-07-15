##don't forget the library dependencies
require(dplyr)
require(ggplot2)
library(ggpubr)
theme_set(theme_pubclean())

msr <- read.csv("/home/daniel/Sync/Google Drive/Niyogi Lab/Zofingiensis/Switch Paper/Switch Paper Part B/Replot/oxygen_consumption.csv", header=TRUE,as.is=TRUE)
#msr <- read.csv("/Users/Meeko/Sync/Google\ Drive/Niyogi\ Lab/Zofingiensis/Switch\ Paper/Switch\ Paper\ Part\ B/Replot/oxygen_consumption.csv", header=TRUE,as.is=TRUE)
msr$strain = factor(msr$strain, levels= c("wt", "hxk1-1", "hxk1-2"))
msr1 <- msr %>% group_by(strain,glc) %>%  dplyr::summarise(average = mean(O2coms), sd = sd(O2coms))
msr$interact <-interaction(msr$strain,msr$glc)
msr2<- right_join(msr,msr1,by=NULL)

p2<- ggplot(msr1, aes(x=strain,y=average,fill=glc,ymin=average-sd,ymax=average+sd))+
  geom_bar(position = "dodge",stat = "identity",color="black",alpha=0.8,width=0.5)+scale_fill_manual(values=c("grey95", "grey15"))+
  geom_point(data=msr,aes(x= strain,y=O2coms,group=msr$interact),inherit.aes = FALSE,shape=21,color="black",fill="white",position = position_dodge(width=0.5),alpha = 0.7,size=4,stroke=1.1)+
  geom_errorbar(width=0.3, position=position_dodge(0.5),size=0.6)

p2
xlabels2 <-c("WT", expression(paste(italic("hxk1-1"))), expression(paste(italic("hxk1-2"))))
p2+scale_x_discrete(labels = xlabels2)+ ylab(label = paste(" O2 Consumption")) + 
  xlab(label = " ")+ theme_classic(base_size = 20)+
  theme(plot.title = element_text(color = "black", size = 30, hjust = 0.5)) +
  scale_y_continuous(breaks = seq(-7,0,1),limits = c(-7,0),expand=c(0,0)) +
  theme(axis.line.x = element_blank(),axis.ticks.x=element_blank())+
  theme(legend.position = "none")
                              
                              
                              