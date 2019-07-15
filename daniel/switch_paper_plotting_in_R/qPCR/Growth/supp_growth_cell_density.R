##don't forget the library dependencies
require(dplyr)
require(ggplot2)
library(ggpubr)
theme_set(theme_classic())



colorlist <-c("darkgreen","orange","grey75","darkgreen","orange","grey75","darkgreen","orange","grey75","darkgreen","orange","grey75","darkgreen","orange","grey75")
colorlist2<-c("darkgreen","blue","grey75")

wt <- read.csv("/Users/Meeko/Sync/Google\ Drive/Niyogi\ Lab/Zofingiensis/Switch\ Paper/Switch\ Paper\ Part\ B/Replot/growth_vol_supp_WT.csv", header=TRUE,as.is=TRUE)
wt1 <- wt %>% group_by(strain,glc,time) %>%  dplyr::summarise(average = mean(vol), sd = sd(vol))
wt2 <-wt1 %>% group_by(time,glc)
wt3<-right_join(wt,wt2,by=NULL)
wt3$interact <-interaction(wt$time,wt$glc)

p <-ggplot(wt3, aes(x=time,y=average, group=interact,color=glc),color=colorlist2) +
  geom_line(size=1)+ scale_color_manual(values=c("darkgreen","blue","grey10"))+
  geom_jitter(aes(x=time,y=vol,color=glc),shape=21,fill="white",width=0.3)+
  geom_errorbar(aes(x=time,ymin=average-sd,ymax=average+sd),color="grey25",alpha = 0.2, width =5)+
  scale_y_continuous(breaks = seq(0,10000000,1000000),limits = c(0,11000000),expand=c(0,0))+
  theme(legend.position=c(0,1), 
        legend.justification=c(0,1),
        legend.direction="vertical")
  



dog6 <- read.csv("/Users/Meeko/Sync/Google\ Drive/Niyogi\ Lab/Zofingiensis/Switch\ Paper/Switch\ Paper\ Part\ B/Replot/growth_vol_supp_hxk1-1.csv", header=TRUE,as.is=TRUE)
dog61 <- dog6 %>% group_by(strain,glc,time) %>%  dplyr::summarise(average = mean(vol), sd = sd(vol))
dog62 <-dog61 %>% group_by(time,glc)
dog63<-right_join(dog6,dog62,by=NULL)
dog63$interact <-interaction(dog6$time,dog6$glc)

p1 <-ggplot(dog63, aes(x=time,y=average, group=interact,color=glc),color=colorlist2) +
  geom_line(size=1)+ scale_color_manual(values=c("darkgreen","blue","grey10"))+
  geom_jitter(aes(x=time,y=vol,color=glc),shape=21,fill="white",width=0.3)+
  geom_errorbar(aes(x=time,ymin=average-sd,ymax=average+sd),color="grey25",alpha = 0.2, width =5)+
  scale_y_continuous(breaks = seq(0,10000000,1000000),limits = c(0,11000000),expand=c(0,0))+
  theme(legend.position="none")


dog9 <- read.csv("/Users/Meeko/Sync/Google\ Drive/Niyogi\ Lab/Zofingiensis/Switch\ Paper/Switch\ Paper\ Part\ B/Replot/growth_vol_supp_hxk1-2.csv", header=TRUE,as.is=TRUE)
dog91 <- dog9 %>% group_by(strain,glc,time) %>%  dplyr::summarise(average = mean(vol), sd = sd(vol))
dog92 <-dog91 %>% group_by(time,glc)
dog93<-right_join(dog9,dog92,by=NULL)
dog93$interact <-interaction(dog9$time,dog9$glc)

p2 <-ggplot(dog93, aes(x=time,y=average, group=interact,color=glc),color=colorlist2) +
  geom_line(size=1)+ scale_color_manual(values=c("darkgreen","blue","grey10"))+
  geom_jitter(aes(x=time,y=vol,color=glc),shape=21,fill="white",width=0.3)+
  geom_errorbar(aes(x=time,ymin=average-sd,ymax=average+sd),color="grey25",alpha = 0.2, width =5)+
  scale_y_continuous(breaks = seq(0,10000000,1000000),limits = c(0,11000000),expand=c(0,0))+
  theme(legend.position="none")

p
ggsave(filename = "wt_Cell_density.pdf", width=3.9, height = 3.273, units="cm", plot=last_plot(),scale = 5 ,device = cairo_pdf,
       path="/Users/Meeko/Sync/Google Drive/Niyogi Lab/Zofingiensis/Switch Paper/Switch Paper Part B/Replot/Plots/To edit/")
p1
ggsave(filename = "hxk1-1_Cell_density.pdf", width=3.9, height = 3.273, units="cm", plot=last_plot(),scale = 5 ,device = cairo_pdf,
       path="/Users/Meeko/Sync/Google Drive/Niyogi Lab/Zofingiensis/Switch Paper/Switch Paper Part B/Replot/Plots/To edit/")
p2
ggsave(filename = "hxk1-2_Cell_density.pdf", width=3.9, height = 3.273, units="cm", plot=last_plot(),scale = 5 ,device = cairo_pdf,
       path="/Users/Meeko/Sync/Google Drive/Niyogi Lab/Zofingiensis/Switch Paper/Switch Paper Part B/Replot/Plots/To edit/")



wt_log <-ggplot(wt3, aes(x=time,y=average, group=interact,color=glc),color=colorlist2) +
  geom_line(size=1)+ scale_color_manual(values=c("darkgreen","blue","grey10"))+
  geom_jitter(aes(x=time,y=vol,color=glc),shape=21,fill="white",width=0.3)+
  geom_errorbar(aes(x=time,ymin=average-sd,ymax=average+sd),color="grey25",alpha = 0.2, width =5)+
  scale_y_log10(limits=c(700000,11000000),expand=c(0,0))+
  theme(legend.position=c(0,1), 
        legend.justification=c(0,1),
        legend.direction="vertical")+
  scale_x_continuous(breaks = c(0,24,48,72,96))

wt_log
  ggsave(filename = "wt_Cell_density_log.pdf", width=3.9, height = 3.273, units="cm", plot=last_plot(),scale = 5 ,device = cairo_pdf,
       path="/Users/Meeko/Sync/Google Drive/Niyogi Lab/Zofingiensis/Switch Paper/Switch Paper Part B/Replot/Plots/To edit/")

d6_log <-ggplot(dog63, aes(x=time,y=average, group=interact,color=glc),color=colorlist2) +
  geom_line(size=1)+ scale_color_manual(values=c("darkgreen","blue","grey10"))+
  geom_jitter(aes(x=time,y=vol,color=glc),shape=21,fill="white",width=0.3)+
  geom_errorbar(aes(x=time,ymin=average-sd,ymax=average+sd),color="grey25",alpha = 0.2, width =5)+
  scale_y_log10(limits=c(700000,11000000),expand=c(0,0))+
  theme(legend.position="none")+
  scale_x_continuous(breaks = c(0,24,48,72,96))
d6_log

ggsave(filename = "hxk1-1_cell_density_log.pdf", width=3.9, height = 3.273, units="cm", plot=last_plot(),scale = 5 ,device = cairo_pdf,
       path="/Users/Meeko/Sync/Google Drive/Niyogi Lab/Zofingiensis/Switch Paper/Switch Paper Part B/Replot/Plots/To edit/")

d9_log <-ggplot(dog93, aes(x=time,y=average, group=interact,color=glc),color=colorlist2) +
  geom_line(size=1)+ scale_color_manual(values=c("darkgreen","blue","grey10"))+
  geom_jitter(aes(x=time,y=vol,color=glc),shape=21,fill="white",width=0.3)+
  geom_errorbar(aes(x=time,ymin=average-sd,ymax=average+sd),color="grey25",alpha = 0.2, width =5)+
  scale_y_log10(limits=c(700000,11000000),expand=c(0,0))+
  theme(legend.position="none")+
  scale_x_continuous(breaks = c(0,24,48,72,96))

d9_log

ggsave(filename = "hxk1-2_cell_density_log.pdf", width=3.9, height = 3.273, units="cm", plot=last_plot(),scale = 5 ,device = cairo_pdf,
       path="/Users/Meeko/Sync/Google Drive/Niyogi Lab/Zofingiensis/Switch Paper/Switch Paper Part B/Replot/Plots/To edit/")
