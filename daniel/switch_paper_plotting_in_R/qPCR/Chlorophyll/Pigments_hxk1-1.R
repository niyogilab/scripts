##don't forget the library dependencies
require(dplyr)
require(ggplot2)
library(ggpubr)
theme_set(theme_pubclean())

#msr <- read.csv("/home/daniel/Sync/Google Drive/Niyogi Lab/Zofingiensis/Switch Paper/Switch Paper Part B/Replot/hxk2_pigments.csv", header=TRUE,as.is=TRUE)
msr <- read.csv("/Users/Meeko/Sync/Google\ Drive/Niyogi\ Lab/Zofingiensis/Switch\ Paper/Switch\ Paper\ Part\ B/Replot/pigments_hxk1-1.csv", header=TRUE,as.is=TRUE)
msr1 <- msr %>% group_by(pigment,glc) %>%  dplyr::summarise(average = mean(nmol), sd = sd(nmol))
msr$interact <-interaction(msr$pigment,msr$glc)
msr2<- right_join(msr,msr1,by=NULL)

p2<- ggplot(msr2, aes(x=pigment,y=average,fill=glc,ymin=average-sd,ymax=average+sd))+
  geom_bar(position = "dodge",stat = "identity",color="black",alpha=0.8,width=0.5)+scale_fill_manual(values=c("grey95", "grey35"))+
  geom_point(data=msr,aes(x= pigment,y=nmol,group=msr$interact),inherit.aes = FALSE,shape=21,color="black",fill="white",position = position_dodge(width=0.5),alpha = 0.7,size=4,stroke=1.1)+
  geom_errorbar(width=0.3, position=position_dodge(0.5),size=0.6)
p2

p2+ 
  xlab(label = " ")+ theme_classic(base_size = 20)+
  theme(plot.title = element_text(color = "black", size = 30, hjust = 0.5)) +
  scale_y_continuous(breaks = seq(0,500,50),limits = c(0,500),expand=c(0,0)) +
  theme(axis.ticks.x=element_blank())+
  theme(legend.position = "none")+
  ylab(NULL)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

ggsave(filename = "Hxk1-1_pigments.pdf", width=5.696, height = 2.765, units="cm", plot=last_plot(),scale = 5 ,path="/Users/Meeko/Sync/Google Drive/Niyogi Lab/Zofingiensis/Switch Paper/Switch Paper Part B/Replot/Plots/To edit/")
