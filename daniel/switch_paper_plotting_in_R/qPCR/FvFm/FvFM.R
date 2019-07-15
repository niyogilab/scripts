##don't forget the library dependencies
require(dplyr)
require(ggplot2)
library(ggpubr)
theme_set(theme_classic())

msr <- read.csv("path/to/FvFm_example_data.csv", header=TRUE,as.is=TRUE)
msr$strain = factor(msr$strain, levels= c("wt", "hxk1-1", "hxk1-2"))
msr1 <- msr %>% group_by(strain,glc) %>%  dplyr::summarise(average = mean(Fv.Fm), sd = sd(Fv.Fm))
msr$interact <-interaction(msr$strain,msr$glc)

#manually jitter points in illustrator
xlabels2 <-c("WT", expression(paste(italic("hxk1"),"-1")), expression(paste(italic("hxk1"),"-2")))
p2<- ggplot(msr1, aes(x=strain,y=average,fill=glc,ymin=average-sd,ymax=average+sd))+
  geom_bar(position = "dodge",stat = "identity",color="black",alpha=0.8,width=0.5)+scale_fill_manual(values=c("grey95", "grey15"))+
  geom_point(data=msr,aes(x= strain,y=Fv.Fm,group=msr$interact),
             inherit.aes = FALSE,shape=21,color="black",
             fill="white",stroke=1.1,position = position_dodge(width=0.5),
             alpha = 0.8, size =4)+
  geom_errorbar(width=0.2, position=position_dodge(0.5),size=0.6)

p2+scale_x_discrete(labels = xlabels2)+ggtitle("Fv/Fm")+ ylab(label = "Fv/Fm") + xlab(label = " ")+ theme_classic(base_size = 20)+
  theme(plot.title = element_text(color = "black", size = 10, hjust = 0.5)) +scale_y_continuous(breaks = seq(0,1,0.1),expand=c(0,0),limits = c(0,0.75))+theme(legend.position = "none")


