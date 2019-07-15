##don't forget the library dependencies
require(dplyr)
require(ggplot2)
library(ggpubr)
theme_set(theme_classic())


dta <- read.csv("/Users/Meeko/Sync/Google\ Drive/Niyogi\ Lab/Zofingiensis/Switch\ Paper/Switch\ Paper\ Part\ B/Replot/supp_groWth_data.csv", header=TRUE,as.is=TRUE)
#weird empty columns

dta$X = NULL
dta$X.1=NULL
dta$X.2=NULL
dta$X.3=NULL
dta$X.4=NULL

colorlist<- c("darkgreen","darkgreen","darkgreen","darkgreen","darkgreen","blue","blue","blue","blue","blue","grey45","grey45","grey45","grey45","grey45")
colorlist2<-c("darkgreen","blue","grey45")

Wt <- dta %>% filter(Strain == "Wt") %>% group_by(Treatment, Hour)
Wt_stat <-Wt %>% group_by(Treatment, Hour)%>%dplyr::summarise(mean.vol.per.ml=mean(vol.per.ml), 
                                                              sd.vol.per.ml = sd(vol.per.ml),
                                                              mean.cell.per.ml=mean(cell.per.ml),
                                                              sd.cell.per.ml = sd(cell.per.ml),
                                                              mean.cell.vol=mean(cell.vol),
                                                              sd.cell.vol=sd(cell.vol))

Wt_plot_cell_ml <- ggplot(data=Wt_stat, aes(x=Hour,y=mean.cell.per.ml, 
                                            ymin = mean.cell.per.ml -sd.cell.per.ml, ymax= mean.cell.per.ml+sd.cell.per.ml))+
  geom_line(stat= "identity",aes(color =Treatment,group= Treatment),size=1)+
  scale_color_manual(values=colorlist2)+
  geom_jitter(data=Wt,aes(x=Hour,y=cell.per.ml,color=Treatment),width = 4,shape=21,fill=NA,inherit.aes = FALSE)+
  geom_errorbar(width = 8, color =colorlist , alpha = 0.8)+
  scale_x_continuous(breaks = c(0,24,48,72,96))+
  theme(legend.position="none")+
  theme(aspect.ratio = 0.8)+
  scale_y_continuous(limits=c(0,10500000),breaks = c(0,1000000,2000000,3000000,4000000,5000000,6000000,7000000,8000000,9000000,10000000),name = NULL,expand=c(0,0))
Wt_plot_cell_ml
ggsave(filename = "wt_cell_per_ml.pdf", plot=last_plot(),width=5.116,height=3.777,units = "cm",scale=3,device = cairo_pdf,
       path="/Users/Meeko/Sync/Google Drive/Niyogi Lab/Zofingiensis/Switch Paper/Switch Paper Part B/Replot/Plots/To edit/Supplemental_growth")

Wt_plot_vol_ml <- ggplot(data=Wt_stat, aes(x=Hour,y=mean.vol.per.ml, 
                                            ymin = mean.vol.per.ml -sd.vol.per.ml, ymax= mean.vol.per.ml+sd.vol.per.ml))+
  geom_line(stat= "identity",aes(color =Treatment,group= Treatment),size=1)+
  scale_color_manual(values=colorlist2)+
  geom_jitter(data=Wt,aes(x=Hour,y=vol.per.ml,color=Treatment),width = 4,shape=21,fill=NA,inherit.aes = FALSE)+
  geom_errorbar(width = 8, color =colorlist , alpha = 0.8)+
  scale_x_continuous(breaks = c(0,24,48,72,96))+
  theme(legend.position="none")+
  theme(aspect.ratio = 0.8)+
  scale_y_continuous(limits=c(0,1233000000),breaks = c(0,2E8,4E8,6E8,8E8,1E9,1.2E9),name = NULL,expand=c(0,0))
Wt_plot_vol_ml

ggsave(filename = "wt_vol_per_ml.pdf", plot=last_plot(),width=5.116,height=3.777,units = "cm",scale=3,device = cairo_pdf,
       path="/Users/Meeko/Sync/Google Drive/Niyogi Lab/Zofingiensis/Switch Paper/Switch Paper Part B/Replot/Plots/To edit/Supplemental_growth")


Wt_plot_cell_vol <- ggplot(data=Wt_stat, aes(x=Hour,y=mean.cell.vol, 
                                           ymin = mean.cell.vol -sd.cell.vol, ymax= mean.cell.vol+sd.cell.vol))+
  geom_line(stat= "identity",aes(color =Treatment,group= Treatment),size=1)+
  scale_color_manual(values=colorlist2)+
  geom_jitter(data=Wt,aes(x=Hour,y=cell.vol,color=Treatment),width = 4,shape=21,fill=NA,inherit.aes = FALSE)+
  geom_errorbar(width = 8, color =colorlist , alpha = 0.8)+
  scale_x_continuous(breaks = c(0,24,48,72,96))+
  theme(legend.position="none")+
  theme(aspect.ratio = 0.8)+
  scale_y_continuous(limits=c(0,300),name = NULL,expand=c(0,0))
Wt_plot_cell_vol

ggsave(filename = "wt_cell_vol.pdf", plot=last_plot(),width=5.116,height=3.777,units = "cm",scale=3,device = cairo_pdf,
       path="/Users/Meeko/Sync/Google Drive/Niyogi Lab/Zofingiensis/Switch Paper/Switch Paper Part B/Replot/Plots/To edit/Supplemental_growth")


################################# -- Hxk1-1 ########################

dog6 <- dta %>% filter(Strain == "Hxk1-1") %>% group_by(Treatment, Hour)
dog6_stat <-dog6 %>% group_by(Treatment, Hour)%>%dplyr::summarise(mean.vol.per.ml=mean(vol.per.ml), 
                                                              sd.vol.per.ml = sd(vol.per.ml),
                                                              mean.cell.per.ml=mean(cell.per.ml),
                                                              sd.cell.per.ml = sd(cell.per.ml),
                                                              mean.cell.vol=mean(cell.vol),
                                                              sd.cell.vol=sd(cell.vol))

dog6_plot_cell_ml <- ggplot(data=dog6_stat, aes(x=Hour,y=mean.cell.per.ml, 
                                            ymin = mean.cell.per.ml -sd.cell.per.ml, ymax= mean.cell.per.ml+sd.cell.per.ml))+
  geom_line(stat= "identity",aes(color =Treatment,group= Treatment),size=1)+
  scale_color_manual(values=colorlist2)+
  geom_jitter(data=dog6,aes(x=Hour,y=cell.per.ml,color=Treatment),width = 4,shape=21,fill=NA,inherit.aes = FALSE)+
  geom_errorbar(width = 8, color =colorlist , alpha = 0.8)+
  scale_x_continuous(breaks = c(0,24,48,72,96))+
  theme(legend.position="none")+
  theme(aspect.ratio = 0.8)+
  scale_y_continuous(limits=c(0,10500000),breaks = c(0,1000000,2000000,3000000,4000000,5000000,6000000,7000000,8000000,9000000,10000000),name = NULL,expand=c(0,0))
dog6_plot_cell_ml
ggsave(filename = "hxk1-1_cell_per_ml.pdf", plot=last_plot(),width=5.116,height=3.777,units = "cm",scale=3,device = cairo_pdf,
       path="/Users/Meeko/Sync/Google Drive/Niyogi Lab/Zofingiensis/Switch Paper/Switch Paper Part B/Replot/Plots/To edit/Supplemental_growth")


dog6_plot_vol_ml <- ggplot(data=dog6_stat, aes(x=Hour,y=mean.vol.per.ml, 
                                           ymin = mean.vol.per.ml -sd.vol.per.ml, ymax= mean.vol.per.ml+sd.vol.per.ml))+
  geom_line(stat= "identity",aes(color =Treatment,group= Treatment),size=1)+
  scale_color_manual(values=colorlist2)+
  geom_jitter(data=dog6,aes(x=Hour,y=vol.per.ml,color=Treatment),width = 4,shape=21,fill=NA,inherit.aes = FALSE)+
  geom_errorbar(width = 8, color =colorlist , alpha = 0.8)+
  scale_x_continuous(breaks = c(0,24,48,72,96))+
  theme(legend.position="none")+
  theme(aspect.ratio = 0.8)+
  scale_y_continuous(limits=c(0,1233000000),breaks = c(0,2E8,4E8,6E8,8E8,1E9,1.2E9),name = NULL,expand=c(0,0))
dog6_plot_vol_ml

ggsave(filename = "hxk1-1_vol_per_ml.pdf", plot=last_plot(),width=5.116,height=3.777,units = "cm",scale=3,device = cairo_pdf,
       path="/Users/Meeko/Sync/Google Drive/Niyogi Lab/Zofingiensis/Switch Paper/Switch Paper Part B/Replot/Plots/To edit/Supplemental_growth")

dog6_plot_cell_vol <- ggplot(data=dog6_stat, aes(x=Hour,y=mean.cell.vol, 
                                             ymin = mean.cell.vol -sd.cell.vol, ymax= mean.cell.vol+sd.cell.vol))+
  geom_line(stat= "identity",aes(color =Treatment,group= Treatment),size=1)+
  scale_color_manual(values=colorlist2)+
  geom_jitter(data=dog6,aes(x=Hour,y=cell.vol,color=Treatment),width = 4,shape=21,fill=NA,inherit.aes = FALSE)+
  geom_errorbar(width = 8, color =colorlist , alpha = 0.8)+
  scale_x_continuous(breaks = c(0,24,48,72,96))+
  theme(legend.position="none")+
  theme(aspect.ratio = 0.8)+
  scale_y_continuous(limits=c(0,300),name = NULL,expand=c(0,0))
dog6_plot_cell_vol

ggsave(filename = "hxk1-1_cell_vol.pdf", plot=last_plot(),width=5.116,height=3.777,units = "cm",scale=3,device = cairo_pdf,
       path="/Users/Meeko/Sync/Google Drive/Niyogi Lab/Zofingiensis/Switch Paper/Switch Paper Part B/Replot/Plots/To edit/Supplemental_growth")


####################################### Hxk1-2

dog9 <- dta %>% filter(Strain == "Hxk1-2") %>% group_by(Treatment, Hour)
dog9_stat <-dog9 %>% group_by(Treatment, Hour)%>%dplyr::summarise(mean.vol.per.ml=mean(vol.per.ml), 
                                                                  sd.vol.per.ml = sd(vol.per.ml),
                                                                  mean.cell.per.ml=mean(cell.per.ml),
                                                                  sd.cell.per.ml = sd(cell.per.ml),
                                                                  mean.cell.vol=mean(cell.vol),
                                                                  sd.cell.vol=sd(cell.vol))

dog9_plot_cell_ml <- ggplot(data=dog9_stat, aes(x=Hour,y=mean.cell.per.ml, 
                                                ymin = mean.cell.per.ml -sd.cell.per.ml, ymax= mean.cell.per.ml+sd.cell.per.ml))+
  geom_line(stat= "identity",aes(color =Treatment,group= Treatment),size=1)+
  scale_color_manual(values=colorlist2)+
  geom_jitter(data=dog9,aes(x=Hour,y=cell.per.ml,color=Treatment),width = 4,shape=21,fill=NA,inherit.aes = FALSE)+
  geom_errorbar(width = 8, color =colorlist , alpha = 0.8)+
  scale_x_continuous(breaks = c(0,24,48,72,96))+
  theme(legend.position="none")+
  theme(aspect.ratio = 0.8)+
  scale_y_continuous(limits=c(0,10500000),breaks = c(0,1000000,2000000,3000000,4000000,5000000,6000000,7000000,8000000,9000000,10000000),name = NULL,expand=c(0,0))
dog9_plot_cell_ml

ggsave(filename = "hxk1-2_cell_per_ml.pdf", plot=last_plot(),width=5.116,height=3.777,units = "cm",scale=3,device = cairo_pdf,
       path="/Users/Meeko/Sync/Google Drive/Niyogi Lab/Zofingiensis/Switch Paper/Switch Paper Part B/Replot/Plots/To edit/Supplemental_growth")


dog9_plot_vol_ml <- ggplot(data=dog9_stat, aes(x=Hour,y=mean.vol.per.ml, 
                                               ymin = mean.vol.per.ml -sd.vol.per.ml, ymax= mean.vol.per.ml+sd.vol.per.ml))+
  geom_line(stat= "identity",aes(color =Treatment,group= Treatment),size=1)+
  scale_color_manual(values=colorlist2)+
  geom_jitter(data=dog9,aes(x=Hour,y=vol.per.ml,color=Treatment),width = 4,shape=21,fill=NA,inherit.aes = FALSE)+
  geom_errorbar(width = 8, color =colorlist , alpha = 0.8)+
  scale_x_continuous(breaks = c(0,24,48,72,96))+
  theme(legend.position="none")+
  theme(aspect.ratio = 0.8)+
  scale_y_continuous(limits=c(0,1233000000),breaks = c(0,2E8,4E8,6E8,8E8,1E9,1.2E9),name = NULL,expand=c(0,0))
dog9_plot_vol_ml

ggsave(filename = "hxk1-2_vol_per_ml.pdf", plot=last_plot(),width=5.116,height=3.777,units = "cm",scale=3,device = cairo_pdf,
       path="/Users/Meeko/Sync/Google Drive/Niyogi Lab/Zofingiensis/Switch Paper/Switch Paper Part B/Replot/Plots/To edit/Supplemental_growth")


dog9_plot_cell_vol <- ggplot(data=dog9_stat, aes(x=Hour,y=mean.cell.vol, 
                                                 ymin = mean.cell.vol -sd.cell.vol, ymax= mean.cell.vol+sd.cell.vol))+
  geom_line(stat= "identity",aes(color =Treatment,group= Treatment),size=1)+
  scale_color_manual(values=colorlist2)+
  geom_jitter(data=dog9,aes(x=Hour,y=cell.vol,color=Treatment),width = 4,shape=21,fill=NA,inherit.aes = FALSE)+
  geom_errorbar(width = 8, color =colorlist , alpha = 0.8)+
  scale_x_continuous(breaks = c(0,24,48,72,96))+
  theme(legend.position="none")+
  theme(aspect.ratio = 0.8)+
  scale_y_continuous(limits=c(0,300),name = NULL,expand=c(0,0))
dog9_plot_cell_vol

ggsave(filename = "hxk1-2_cell_vol.pdf", plot=last_plot(),width=5.116,height=3.777,units = "cm",scale=3,device = cairo_pdf,
       path="/Users/Meeko/Sync/Google Drive/Niyogi Lab/Zofingiensis/Switch Paper/Switch Paper Part B/Replot/Plots/To edit/Supplemental_growth")

dev.off()

