##don't forget the library dependencies
require(dplyr)
require(ggplot2)
library(ggpubr)
theme_set(theme_pubclean())

#read in the data. 

qpcr <- read.csv("/home/daniel//Sync/Google Drive/Niyogi Lab/Zofingiensis/Switch Paper/qPCR/qPCR_Data_all2.csv", header=TRUE,as.is=TRUE) %>% tbl_df
#average the technnical replicates
qpcr <- qpcr %>% group_by(strain, bio_reps, treatment, timepoint, rxn, target) %>% mutate(mean_tech_rep_ct = mean(ct)) %>% ungroup
#remove redundant data. The tecnical reps are averaged, so there isn't any reason to keep 3 of them moving forward. 
qpcr2 <-qpcr %>% filter(tech_reps ==1)
#all work will now be on the averaged tech reps, so raw CT can be removed. 
qpcr2$ct = NULL
qpcr2$tech_reps = NULL

#APPS is the control gene to calculate Delta CT. 
qpcr_apps <- qpcr2 %>% filter(target == "APPS") %>% group_by(strain, bio_reps, treatment, timepoint, rxn) %>% summarize(apps_ct = mean(mean_tech_rep_ct))
qpcr3 <- left_join(qpcr2, qpcr_apps, by=c("strain", "bio_reps", "treatment", "timepoint", "rxn")) %>% filter(target != "APPS")
qpcr4 <- qpcr3 %>% mutate(delta_ct = mean_tech_rep_ct - apps_ct)

DeltaCT <- qpcr4 %>% distinct(strain, treatment, target, timepoint, bio_reps, delta_ct)

#reorder the data to have a treatment and control CT column

control<- DeltaCT %>% filter(treatment =="(-)")
control <-dplyr::rename(control, control_ct = delta_ct)

treatment <-DeltaCT %>% filter(treatment =="(+)")
treatment <-dplyr::rename(treatment, treatment_ct = delta_ct)

control$treatment_ct = treatment$treatment_ct 
reorder <-control
reorder$treatment = NULL

#calculate delta delta CT = treatment-control
FoldChange <- mutate(reorder, ddCT = treatment_ct-control_ct)

#log2 FC


FoldChange$fc = 2^-(FoldChange$ddCT)
FoldChange$lg2FC = log2(FoldChange$fc)
FoldChange_stat <- FoldChange %>% group_by(strain,target,timepoint)%>% dplyr::summarise(mean_lg2_FC=mean(lg2FC), sd = sd(lg2FC))

##plotting 0.5h

xlabels <-c("WT",expression(paste(italic("hxk1-1"))), expression(paste(italic("hxk1-2"))))

hxk0.5<-ggplot(data=FoldChange_stat %>% filter(target =="HXK1", timepoint ==0.5), 
               aes(x=strain, y=mean_lg2_FC,ymin= mean_lg2_FC-sd, ymax = mean_lg2_FC+sd))+
  geom_bar(stat="identity",position = "dodge",width=0.3, fill="grey95",color="black")+
  geom_errorbar(stat="identity",width = 0.15, color = "black")+
  geom_jitter(data=FoldChange %>%filter(target =="HXK1", timepoint ==0.5), 
              aes(x=strain,y=lg2FC),inherit.aes = FALSE, width = 0.08, size = 3, shape = 21, fill = "white")+
  scale_x_discrete(labels = xlabels)+
  scale_y_continuous(limits = c(-8,10), expand=c(0,0),breaks = seq(-8,10,2))+
  geom_hline(yintercept = 0,size=0.5) +
  ylab(bquote('Log'[2]~ 'Fold Change at 0.5 h'))+
  xlab(NULL)+
  theme_classic()+
  ggtitle(expression(paste(italic("HXK1"))))+
  theme(plot.title = element_text(hjust = 0.5,size =25))+
  theme(axis.title.y = element_text(size = 20),axis.text=element_text(size=16),axis.line.x = element_blank(),axis.ticks.x=element_blank())+
  theme(axis.text.x = element_text(angle = 40, hjust = 1))

hxk0.5
ggsave(filename = "HXK1_qPCR_0.5.pdf", plot=last_plot(),width=7,height=8.5,units = "cm",scale=2,useDingbats=FALSE,
       path="/Users/Meeko/Sync/Google Drive/Niyogi Lab/Zofingiensis/Switch Paper/Switch Paper Part B/Replot/Plots/To edit/qPCR")

FAD20.5<-ggplot(data=FoldChange_stat %>% filter(target =="FAD2", timepoint ==0.5), 
                aes(x=strain, y=mean_lg2_FC,ymin= mean_lg2_FC-sd, ymax = mean_lg2_FC+sd))+
  geom_bar(stat="identity",position = "dodge",width=0.3, fill="grey95",color="black")+
  geom_errorbar(stat="identity",width = 0.15, color = "black")+
  geom_jitter(data=FoldChange %>%filter(target =="FAD2", timepoint ==0.5), 
              aes(x=strain,y=lg2FC),inherit.aes = FALSE, width = 0.08, size = 3, shape = 21, fill = "white")+
  scale_x_discrete(labels = xlabels)+
  scale_y_continuous(limits = c(-8,10), expand=c(0,0),breaks = seq(-8,10,2))+
  geom_hline(yintercept = 0,size=0.5) +
  ylab(bquote('Log'[2]~ 'Fold Change at 0.5 h'))+
  xlab(NULL)+
  theme_classic()+
  ggtitle(expression(paste(italic("FAD2"))))+
  theme(plot.title = element_text(hjust = 0.5,size =25))+
  theme(axis.title.y = element_text(size = 20),axis.text=element_text(size=16),axis.line.x = element_blank(),axis.ticks.x=element_blank())+
  theme(axis.text.x = element_text(angle = 40, hjust = 1))

FAD20.5
  ggsave(filename = "FAD2_qPCR_0.5.pdf", plot=last_plot(),width=7,height=8.5,units = "cm",scale=2,
       path="/Users/Meeko/Sync/Google Drive/Niyogi Lab/Zofingiensis/Switch Paper/Switch Paper Part B/Replot/Plots/To edit/qPCR")


GLK10.5<-ggplot(data=FoldChange_stat %>% filter(target =="GLK1", timepoint ==0.5), 
                aes(x=strain, y=mean_lg2_FC,ymin= mean_lg2_FC-sd, ymax = mean_lg2_FC+sd))+
  geom_bar(stat="identity",position = "dodge",width=0.3, fill="grey95",color="black")+
  geom_errorbar(stat="identity",width = 0.15, color = "black")+
  geom_jitter(data=FoldChange %>%filter(target =="GLK1", timepoint ==0.5), 
              aes(x=strain,y=lg2FC),inherit.aes = FALSE, width = 0.08, size = 3, shape = 21, fill = "white")+
  scale_x_discrete(labels = xlabels)+
  scale_y_continuous(limits = c(-8,10), expand=c(0,0),breaks = seq(-8,10,2))+
  geom_hline(yintercept = 0,size=0.5) +
  ylab(bquote('Log'[2]~ 'Fold Change at 0.5 h'))+
  xlab(NULL)+
  theme_classic()+
  ggtitle(expression(paste(italic("GLK1"))))+
  theme(plot.title = element_text(hjust = 0.5,size =25))+
  theme(axis.title.y = element_text(size = 20),axis.text=element_text(size=16),axis.line.x = element_blank(),axis.ticks.x=element_blank())+
  theme(axis.text.x = element_text(angle = 40, hjust = 1))

GLK10.5
ggsave(filename = "GLK1_qPCR_0.5.pdf", plot=last_plot(),width=7,height=8.5,units = "cm",scale=2,
       path="/Users/Meeko/Sync/Google Drive/Niyogi Lab/Zofingiensis/Switch Paper/Switch Paper Part B/Replot/Plots/To edit/qPCR")

GAP10.5<-ggplot(data=FoldChange_stat %>% filter(target =="GAP1", timepoint ==0.5), 
                aes(x=strain, y=mean_lg2_FC,ymin= mean_lg2_FC-sd, ymax = mean_lg2_FC+sd))+
  geom_bar(stat="identity",position = "dodge",width=0.3, fill="grey95",color="black")+
  geom_errorbar(stat="identity",width = 0.15, color = "black")+
  geom_jitter(data=FoldChange %>%filter(target =="GAP1", timepoint ==0.5), 
              aes(x=strain,y=lg2FC),inherit.aes = FALSE, width = 0.08, size = 3, shape = 21, fill = "white")+
  scale_x_discrete(labels = xlabels)+
  scale_y_continuous(limits = c(-8,10), expand=c(0,0),breaks = seq(-8,10,2))+
  geom_hline(yintercept = 0,size=0.5) +
  ylab(bquote('Log'[2]~ 'Fold Change at 0.5 h'))+
  xlab(NULL)+
  theme_classic()+
  ggtitle(expression(paste(italic("GAP1"))))+
  theme(plot.title = element_text(hjust = 0.5,size =25))+
  theme(axis.title.y = element_text(size = 20),axis.text=element_text(size=16),axis.line.x = element_blank(),axis.ticks.x=element_blank())+
  theme(axis.text.x = element_text(angle = 40, hjust = 1))

GAP10.5
ggsave(filename = "GAP1_qPCR_0.5.pdf", plot=last_plot(),width=7,height=8.5,units = "cm",scale=2,
       path="/Users/Meeko/Sync/Google Drive/Niyogi Lab/Zofingiensis/Switch Paper/Switch Paper Part B/Replot/Plots/To edit/qPCR")

BKT10.5<-ggplot(data=FoldChange_stat %>% filter(target =="BKT1", timepoint ==0.5), 
                aes(x=strain, y=mean_lg2_FC,ymin= mean_lg2_FC-sd, ymax = mean_lg2_FC+sd))+
  geom_bar(stat="identity",position = "dodge",width=0.3, fill="grey95",color="black")+
  geom_errorbar(stat="identity",width = 0.15, color = "black")+
  geom_jitter(data=FoldChange %>%filter(target =="BKT1", timepoint ==0.5), 
              aes(x=strain,y=lg2FC),inherit.aes = FALSE, width = 0.08, size = 3, shape = 21, fill = "white")+
  scale_x_discrete(labels = xlabels)+
  scale_y_continuous(limits = c(-8,10), expand=c(0,0),breaks = seq(-8,10,2))+
  geom_hline(yintercept = 0,size=0.5) +
  ylab(bquote('Log'[2]~ 'Fold Change at 0.5 h'))+
  xlab(NULL)+
  theme_classic()+
  ggtitle(expression(paste(italic("BKT1"))))+
  theme(plot.title = element_text(hjust = 0.5,size =25))+
  theme(axis.title.y = element_text(size = 20),axis.text=element_text(size=16),axis.line.x = element_blank(),axis.ticks.x=element_blank())+
  theme(axis.text.x = element_text(angle = 40, hjust = 1))

BKT10.5
ggsave(filename = "BKT1_qPCR_0.5.pdf", plot=last_plot(),width=7,height=8.5,units = "cm",scale=2,
       path="/Users/Meeko/Sync/Google Drive/Niyogi Lab/Zofingiensis/Switch Paper/Switch Paper Part B/Replot/Plots/To edit/qPCR")
PSAH0.5<-ggplot(data=FoldChange_stat %>% filter(target =="PSAH", timepoint ==0.5), 
                aes(x=strain, y=mean_lg2_FC,ymin= mean_lg2_FC-sd, ymax = mean_lg2_FC+sd))+
  geom_bar(stat="identity",position = "dodge",width=0.3, fill="grey95",color="black")+
  geom_errorbar(stat="identity",width = 0.15, color = "black")+
  geom_jitter(data=FoldChange %>%filter(target =="PSAH", timepoint ==0.5), 
              aes(x=strain,y=lg2FC),inherit.aes = FALSE, width = 0.08, size = 3, shape = 21, fill = "white")+
  scale_x_discrete(labels = xlabels)+
  scale_y_continuous(limits = c(-8,10), expand=c(0,0),breaks = seq(-8,10,2))+
  geom_hline(yintercept = 0,size=0.5) +
  ylab(bquote('Log'[2]~ 'Fold Change at 0.5 h'))+
  xlab(NULL)+
  theme_classic()+
  ggtitle(expression(paste(italic("PSAH"))))+
  theme(plot.title = element_text(hjust = 0.5,size =25))+
  theme(axis.title.y = element_text(size = 20),axis.text=element_text(size=16),axis.line.x = element_blank(),axis.ticks.x=element_blank())+
  theme(axis.text.x = element_text(angle = 40, hjust = 1))

PSAH0.5
ggsave(filename = "PSAH_qPCR_0.5.pdf", plot=last_plot(),width=7,height=8.5,units = "cm",scale=2,
       path="/Users/Meeko/Sync/Google Drive/Niyogi Lab/Zofingiensis/Switch Paper/Switch Paper Part B/Replot/Plots/To edit/qPCR")

PsbO0.5<-ggplot(data=FoldChange_stat %>% filter(target =="PsbO", timepoint ==0.5), 
                aes(x=strain, y=mean_lg2_FC,ymin= mean_lg2_FC-sd, ymax = mean_lg2_FC+sd))+
  geom_bar(stat="identity",position = "dodge",width=0.3, fill="grey95",color="black")+
  geom_errorbar(stat="identity",width = 0.15, color = "black")+
  geom_jitter(data=FoldChange %>%filter(target =="PsbO", timepoint ==0.5), 
              aes(x=strain,y=lg2FC),inherit.aes = FALSE, width = 0.08, size = 3, shape = 21, fill = "white")+
  scale_x_discrete(labels = xlabels)+
  scale_y_continuous(limits = c(-8,10), expand=c(0,0),breaks = seq(-8,10,2))+
  geom_hline(yintercept = 0,size=0.5) +
  ylab(bquote('Log'[2]~ 'Fold Change at 0.5 h'))+
  xlab(NULL)+
  theme_classic()+
  ggtitle(expression(paste(italic("PsbO"))))+
  theme(plot.title = element_text(hjust = 0.5,size =25))+
  theme(axis.title.y = element_text(size = 20),axis.text=element_text(size=16),axis.line.x = element_blank(),axis.ticks.x=element_blank())+
  theme(axis.text.x = element_text(angle = 40, hjust = 1))

PsbO0.5
ggsave(filename = "PsbO_qPCR_0.5.pdf", plot=last_plot(),width=7,height=8.5,units = "cm",scale=2,
       path="/Users/Meeko/Sync/Google Drive/Niyogi Lab/Zofingiensis/Switch Paper/Switch Paper Part B/Replot/Plots/To edit/qPCR")

RBCS10.5<-ggplot(data=FoldChange_stat %>% filter(target =="RBCS1", timepoint ==0.5), 
                 aes(x=strain, y=mean_lg2_FC,ymin= mean_lg2_FC-sd, ymax = mean_lg2_FC+sd))+
  geom_bar(stat="identity",position = "dodge",width=0.3, fill="grey95",color="black")+
  geom_errorbar(stat="identity",width = 0.15, color = "black")+
  geom_jitter(data=FoldChange %>%filter(target =="RBCS1", timepoint ==0.5), 
              aes(x=strain,y=lg2FC),inherit.aes = FALSE, width = 0.08, size = 3, shape = 21, fill = "white")+
  scale_x_discrete(labels = xlabels)+
  scale_y_continuous(limits = c(-8,10), expand=c(0,0),breaks = seq(-8,10,2))+
  geom_hline(yintercept = 0,size=0.5) +
  ylab(bquote('Log'[2]~ 'Fold Change at 0.5 h'))+
  xlab(NULL)+
  theme_classic()+
  ggtitle(expression(paste(italic("RBCS1"))))+
  theme(plot.title = element_text(hjust = 0.5,size =25))+
  theme(axis.title.y = element_text(size = 20),axis.text=element_text(size=16),axis.line.x = element_blank(),axis.ticks.x=element_blank())+
  theme(axis.text.x = element_text(angle = 40, hjust = 1))

RBCS10.5
ggsave(filename = "RBCS1_qPCR_0.5.pdf", plot=last_plot(),width=7,height=8.5,units = "cm",scale=2,
       path="/Users/Meeko/Sync/Google Drive/Niyogi Lab/Zofingiensis/Switch Paper/Switch Paper Part B/Replot/Plots/To edit/qPCR")

MLDP10.5<-ggplot(data=FoldChange_stat %>% filter(target =="MLDP1", timepoint ==0.5), 
                 aes(x=strain, y=mean_lg2_FC,ymin= mean_lg2_FC-sd, ymax = mean_lg2_FC+sd))+
  geom_bar(stat="identity",position = "dodge",width=0.3, fill="grey95",color="black")+
  geom_errorbar(stat="identity",width = 0.15, color = "black")+
  geom_jitter(data=FoldChange %>%filter(target =="MLDP1", timepoint ==0.5), 
              aes(x=strain,y=lg2FC),inherit.aes = FALSE, width = 0.08, size = 3, shape = 21, fill = "white")+
  scale_x_discrete(labels = xlabels)+
  scale_y_continuous(limits = c(-8,10), expand=c(0,0),breaks = seq(-8,10,2))+
  geom_hline(yintercept = 0,size=0.5) +
  ylab(bquote('Log'[2]~ 'Fold Change at 0.5 h'))+
  xlab(NULL)+
  theme_classic()+
  ggtitle(expression(paste(italic("MLDP1"))))+
  theme(plot.title = element_text(hjust = 0.5,size =25))+
  theme(axis.title.y = element_text(size = 20),axis.text=element_text(size=16),axis.line.x = element_blank(),axis.ticks.x=element_blank())+
  theme(axis.text.x = element_text(angle = 40, hjust = 1))

MLDP10.5
ggsave(filename = "MLDP1_qPCR_0.5.pdf", plot=last_plot(),width=7,height=8.5,units = "cm",scale=2,
       path="/Users/Meeko/Sync/Google Drive/Niyogi Lab/Zofingiensis/Switch Paper/Switch Paper Part B/Replot/Plots/To edit/qPCR")

SMP10.5<-ggplot(data=FoldChange_stat %>% filter(target =="SMP1", timepoint ==0.5), 
                aes(x=strain, y=mean_lg2_FC,ymin= mean_lg2_FC-sd, ymax = mean_lg2_FC+sd))+
  geom_bar(stat="identity",position = "dodge",width=0.3, fill="grey95",color="black")+
  geom_errorbar(stat="identity",width = 0.15, color = "black")+
  geom_jitter(data=FoldChange %>%filter(target =="SMP1", timepoint ==0.5), 
              aes(x=strain,y=lg2FC),inherit.aes = FALSE, width = 0.08, size = 3, shape = 21, fill = "white")+
  scale_x_discrete(labels = xlabels)+
  scale_y_continuous(limits = c(-8,10), expand=c(0,0),breaks = seq(-8,10,2))+
  geom_hline(yintercept = 0,size=0.5) +
  ylab(bquote('Log'[2]~ 'Fold Change at 0.5 h'))+
  xlab(NULL)+
  theme_classic()+
  ggtitle(expression(paste(italic("SMP1"))))+
  theme(plot.title = element_text(hjust = 0.5,size =25))+
  theme(axis.title.y = element_text(size = 20),axis.text=element_text(size=16),axis.line.x = element_blank(),axis.ticks.x=element_blank())+
  theme(axis.text.x = element_text(angle = 40, hjust = 1))

SMP10.5
ggsave(filename = "SMP1_qPCR_0.5.pdf", plot=last_plot(),width=7,height=8.5,units = "cm",scale=2,
       path="/Users/Meeko/Sync/Google Drive/Niyogi Lab/Zofingiensis/Switch Paper/Switch Paper Part B/Replot/Plots/To edit/qPCR")

PYK20.5<-ggplot(data=FoldChange_stat %>% filter(target =="PYK2", timepoint ==0.5), 
                aes(x=strain, y=mean_lg2_FC,ymin= mean_lg2_FC-sd, ymax = mean_lg2_FC+sd))+
  geom_bar(stat="identity",position = "dodge",width=0.3, fill="grey95",color="black")+
  geom_errorbar(stat="identity",width = 0.15, color = "black")+
  geom_jitter(data=FoldChange %>%filter(target =="PYK2", timepoint ==0.5), 
              aes(x=strain,y=lg2FC),inherit.aes = FALSE, width = 0.08, size = 3, shape = 21, fill = "white")+
  scale_x_discrete(labels = xlabels)+
  scale_y_continuous(limits = c(-8,10), expand=c(0,0),breaks = seq(-8,10,2))+
  geom_hline(yintercept = 0,size=0.5) +
  ylab(bquote('Log'[2]~ 'Fold Change at 0.5 h'))+
  xlab(NULL)+
  theme_classic()+
  ggtitle(expression(paste(italic("PYK2"))))+
  theme(plot.title = element_text(hjust = 0.5,size =25))+
  theme(axis.title.y = element_text(size = 20),axis.text=element_text(size=16),axis.line.x = element_blank(),axis.ticks.x=element_blank())+
  theme(axis.text.x = element_text(angle = 40, hjust = 1))

PYK20.5
ggsave(filename = "PYK2_qPCR_0.5.pdf", plot=last_plot(),width=7,height=8.5,units = "cm",scale=2,
       path="/Users/Meeko/Sync/Google Drive/Niyogi Lab/Zofingiensis/Switch Paper/Switch Paper Part B/Replot/Plots/To edit/qPCR")

LHC160.5<-ggplot(data=FoldChange_stat %>% filter(target =="LHC16", timepoint ==0.5), 
                 aes(x=strain, y=mean_lg2_FC,ymin= mean_lg2_FC-sd, ymax = mean_lg2_FC+sd))+
  geom_bar(stat="identity",position = "dodge",width=0.3, fill="grey95",color="black")+
  geom_errorbar(stat="identity",width = 0.15, color = "black")+
  geom_jitter(data=FoldChange %>%filter(target =="LHC16", timepoint ==0.5), 
              aes(x=strain,y=lg2FC),inherit.aes = FALSE, width = 0.08, size = 3, shape = 21, fill = "white")+
  scale_x_discrete(labels = xlabels)+
  scale_y_continuous(limits = c(-8,10), expand=c(0,0),breaks = seq(-8,10,2))+
  geom_hline(yintercept = 0,size=0.5) +
  ylab(bquote('Log'[2]~ 'Fold Change at 0.5 h'))+
  xlab(NULL)+
  theme_classic()+
  ggtitle(expression(paste(italic("LHC16"))))+
  theme(plot.title = element_text(hjust = 0.5,size =25))+
  theme(axis.title.y = element_text(size = 20),axis.text=element_text(size=16),axis.line.x = element_blank(),axis.ticks.x=element_blank())+
  theme(axis.text.x = element_text(angle = 40, hjust = 1))

LHC160.5
ggsave(filename = "LHC16_qPCR_0.5.pdf", plot=last_plot(),width=7,height=8.5,units = "cm",scale=2,
       path="/Users/Meeko/Sync/Google Drive/Niyogi Lab/Zofingiensis/Switch Paper/Switch Paper Part B/Replot/Plots/To edit/qPCR")



COX100.5<-ggplot(data=FoldChange_stat %>% filter(target =="COX10", timepoint ==0.5), 
                 aes(x=strain, y=mean_lg2_FC,ymin= mean_lg2_FC-sd, ymax = mean_lg2_FC+sd))+
  geom_bar(stat="identity",position = "dodge",width=0.3, fill="grey95",color="black")+
  geom_errorbar(stat="identity",width = 0.15, color = "black")+
  geom_jitter(data=FoldChange %>%filter(target =="COX10", timepoint ==0.5), 
              aes(x=strain,y=lg2FC),inherit.aes = FALSE, width = 0.08, size = 3, shape = 21, fill = "white")+
  scale_x_discrete(labels = xlabels)+
  scale_y_continuous(limits = c(-8,10), expand=c(0,0),breaks = seq(-8,10,2))+
  geom_hline(yintercept = 0,size=0.5) +
  ylab(bquote('Log'[2]~ 'Fold Change at 0.5 h'))+
  xlab(NULL)+
  theme_classic()+
  ggtitle(expression(paste(italic("COX10"))))+
  theme(plot.title = element_text(hjust = 0.5,size =25))+
  theme(axis.title.y = element_text(size = 20),axis.text=element_text(size=16),axis.line.x = element_blank(),axis.ticks.x=element_blank())+
  theme(axis.text.x = element_text(angle = 40, hjust = 1))

COX100.5
ggsave(filename = "COX10_qPCR_0.5.pdf", plot=last_plot(),width=7,height=8.5,units = "cm",scale=2,
       path="/Users/Meeko/Sync/Google Drive/Niyogi Lab/Zofingiensis/Switch Paper/Switch Paper Part B/Replot/Plots/To edit/qPCR")
SAD10.5<-ggplot(data=FoldChange_stat %>% filter(target =="SAD1", timepoint ==0.5), 
                aes(x=strain, y=mean_lg2_FC,ymin= mean_lg2_FC-sd, ymax = mean_lg2_FC+sd))+
  geom_bar(stat="identity",position = "dodge",width=0.3, fill="grey95",color="black")+
  geom_errorbar(stat="identity",width = 0.15, color = "black")+
  geom_jitter(data=FoldChange %>%filter(target =="SAD1", timepoint ==0.5), 
              aes(x=strain,y=lg2FC),inherit.aes = FALSE, width = 0.08, size = 3, shape = 21, fill = "white")+
  scale_x_discrete(labels = xlabels)+
  scale_y_continuous(limits = c(-8,10), expand=c(0,0),breaks = seq(-8,10,2))+
  geom_hline(yintercept = 0,size=0.5) +
  ylab(bquote('Log'[2]~ 'Fold Change at 0.5 h'))+
  xlab(NULL)+
  theme_classic()+
  ggtitle(expression(paste(italic("SAD1"))))+
  theme(plot.title = element_text(hjust = 0.5,size =25))+
  theme(axis.title.y = element_text(size = 20),axis.text=element_text(size=16),axis.line.x = element_blank(),axis.ticks.x=element_blank())+
  theme(axis.text.x = element_text(angle = 40, hjust = 1))

SAD10.5
ggsave(filename = "SAD1_qPCR_0.5.pdf", plot=last_plot(),width=7,height=8.5,units = "cm",scale=2,
       path="/Users/Meeko/Sync/Google Drive/Niyogi Lab/Zofingiensis/Switch Paper/Switch Paper Part B/Replot/Plots/To edit/qPCR")

##plotting 12h

xlabels <-c("WT",expression(paste(italic("hxk1-1"))), expression(paste(italic("hxk1-2"))))

hxk12<-ggplot(data=FoldChange_stat %>% filter(target =="HXK1", timepoint ==12), 
              aes(x=strain, y=mean_lg2_FC,ymin= mean_lg2_FC-sd, ymax = mean_lg2_FC+sd))+
  geom_bar(stat="identity",position = "dodge",width=0.3, fill="grey35",color="black")+
  geom_errorbar(stat="identity",width = 0.15, color = "black")+
  geom_jitter(data=FoldChange %>%filter(target =="HXK1", timepoint ==12), 
              aes(x=strain,y=lg2FC),inherit.aes = FALSE, width = 0.08, size = 3, shape = 21, fill = "white")+
  scale_x_discrete(labels = xlabels)+
  scale_y_continuous(limits = c(-8,10), expand=c(0,0),breaks = seq(-8,10,2))+
  geom_hline(yintercept = 0,size=0.5) +
  ylab(bquote('Log'[2]~ 'Fold Change at 12 h'))+
  xlab(NULL)+
  theme_classic()+
  ggtitle(expression(paste(italic("HXK1"))))+
  theme(plot.title = element_text(hjust = 0.5,size =25))+
  theme(axis.title.y = element_text(size = 20),axis.text=element_text(size=16),axis.line.x = element_blank(),axis.ticks.x=element_blank())+
  theme(axis.text.x = element_text(angle = 40, hjust = 1))

hxk12
ggsave(filename = "HXK1_qPCR_12.pdf", plot=last_plot(),width=7,height=8.5,units = "cm",scale=2,
       path="/Users/Meeko/Sync/Google Drive/Niyogi Lab/Zofingiensis/Switch Paper/Switch Paper Part B/Replot/Plots/To edit/qPCR")

FAD212<-ggplot(data=FoldChange_stat %>% filter(target =="FAD2", timepoint ==12), 
               aes(x=strain, y=mean_lg2_FC,ymin= mean_lg2_FC-sd, ymax = mean_lg2_FC+sd))+
  geom_bar(stat="identity",position = "dodge",width=0.3, fill="grey35",color="black")+
  geom_errorbar(stat="identity",width = 0.15, color = "black")+
  geom_jitter(data=FoldChange %>%filter(target =="FAD2", timepoint ==12), 
              aes(x=strain,y=lg2FC),inherit.aes = FALSE, width = 0.08, size = 3, shape = 21, fill = "white")+
  scale_x_discrete(labels = xlabels)+
  scale_y_continuous(limits = c(-8,10), expand=c(0,0),breaks = seq(-8,10,2))+
  geom_hline(yintercept = 0,size=0.5) +
  ylab(bquote('Log'[2]~ 'Fold Change at 12 h'))+
  xlab(NULL)+
  theme_classic()+
  ggtitle(expression(paste(italic("FAD2"))))+
  theme(plot.title = element_text(hjust = 0.5,size =25))+
  theme(axis.title.y = element_text(size = 20),axis.text=element_text(size=16),axis.line.x = element_blank(),axis.ticks.x=element_blank())+
  theme(axis.text.x = element_text(angle = 40, hjust = 1))

FAD212
ggsave(filename = "FAD2_qPCR_12.pdf", plot=last_plot(),width=7,height=8.5,units = "cm",scale=2,
       path="/Users/Meeko/Sync/Google Drive/Niyogi Lab/Zofingiensis/Switch Paper/Switch Paper Part B/Replot/Plots/To edit/qPCR")


GLK112<-ggplot(data=FoldChange_stat %>% filter(target =="GLK1", timepoint ==12), 
               aes(x=strain, y=mean_lg2_FC,ymin= mean_lg2_FC-sd, ymax = mean_lg2_FC+sd))+
  geom_bar(stat="identity",position = "dodge",width=0.3, fill="grey35",color="black")+
  geom_errorbar(stat="identity",width = 0.15, color = "black")+
  geom_jitter(data=FoldChange %>%filter(target =="GLK1", timepoint ==12), 
              aes(x=strain,y=lg2FC),inherit.aes = FALSE, width = 0.08, size = 3, shape = 21, fill = "white")+
  scale_x_discrete(labels = xlabels)+
  scale_y_continuous(limits = c(-8,10), expand=c(0,0),breaks = seq(-8,10,2))+
  geom_hline(yintercept = 0,size=0.5) +
  ylab(bquote('Log'[2]~ 'Fold Change at 12 h'))+
  xlab(NULL)+
  theme_classic()+
  ggtitle(expression(paste(italic("GLK1"))))+
  theme(plot.title = element_text(hjust = 0.5,size =25))+
  theme(axis.title.y = element_text(size = 20),axis.text=element_text(size=16),axis.line.x = element_blank(),axis.ticks.x=element_blank())+
  theme(axis.text.x = element_text(angle = 40, hjust = 1))

GLK112
ggsave(filename = "GLK1_qPCR_12.pdf", plot=last_plot(),width=7,height=8.5,units = "cm",scale=2,
       path="/Users/Meeko/Sync/Google Drive/Niyogi Lab/Zofingiensis/Switch Paper/Switch Paper Part B/Replot/Plots/To edit/qPCR")

GAP112<-ggplot(data=FoldChange_stat %>% filter(target =="GAP1", timepoint ==12), 
               aes(x=strain, y=mean_lg2_FC,ymin= mean_lg2_FC-sd, ymax = mean_lg2_FC+sd))+
  geom_bar(stat="identity",position = "dodge",width=0.3, fill="grey35",color="black")+
  geom_errorbar(stat="identity",width = 0.15, color = "black")+
  geom_jitter(data=FoldChange %>%filter(target =="GAP1", timepoint ==12), 
              aes(x=strain,y=lg2FC),inherit.aes = FALSE, width = 0.08, size = 3, shape = 21, fill = "white")+
  scale_x_discrete(labels = xlabels)+
  scale_y_continuous(limits = c(-8,10), expand=c(0,0),breaks = seq(-8,10,2))+
  geom_hline(yintercept = 0,size=0.5) +
  ylab(bquote('Log'[2]~ 'Fold Change at 12 h'))+
  xlab(NULL)+
  theme_classic()+
  ggtitle(expression(paste(italic("GAP1"))))+
  theme(plot.title = element_text(hjust = 0.5,size =25))+
  theme(axis.title.y = element_text(size = 20),axis.text=element_text(size=16),axis.line.x = element_blank(),axis.ticks.x=element_blank())+
  theme(axis.text.x = element_text(angle = 40, hjust = 1))

GAP112
ggsave(filename = "GAP1_qPCR_12.pdf", plot=last_plot(),width=7,height=8.5,units = "cm",scale=2,
       path="/Users/Meeko/Sync/Google Drive/Niyogi Lab/Zofingiensis/Switch Paper/Switch Paper Part B/Replot/Plots/To edit/qPCR")

BKT112<-ggplot(data=FoldChange_stat %>% filter(target =="BKT1", timepoint ==12), 
               aes(x=strain, y=mean_lg2_FC,ymin= mean_lg2_FC-sd, ymax = mean_lg2_FC+sd))+
  geom_bar(stat="identity",position = "dodge",width=0.3, fill="grey35",color="black")+
  geom_errorbar(stat="identity",width = 0.15, color = "black")+
  geom_jitter(data=FoldChange %>%filter(target =="BKT1", timepoint ==12), 
              aes(x=strain,y=lg2FC),inherit.aes = FALSE, width = 0.08, size = 3, shape = 21, fill = "white")+
  scale_x_discrete(labels = xlabels)+
  scale_y_continuous(limits = c(-8,10), expand=c(0,0),breaks = seq(-8,10,2))+
  geom_hline(yintercept = 0,size=0.5) +
  ylab(bquote('Log'[2]~ 'Fold Change at 12 h'))+
  xlab(NULL)+
  theme_classic()+
  ggtitle(expression(paste(italic("BKT1"))))+
  theme(plot.title = element_text(hjust = 0.5,size =25))+
  theme(axis.title.y = element_text(size = 20),axis.text=element_text(size=16),axis.line.x = element_blank(),axis.ticks.x=element_blank())+
  theme(axis.text.x = element_text(angle = 40, hjust = 1))
BKT112
ggsave(filename = "BKT1_qPCR_12.pdf", plot=last_plot(),width=7,height=8.5,units = "cm",scale=2,
       path="/Users/Meeko/Sync/Google Drive/Niyogi Lab/Zofingiensis/Switch Paper/Switch Paper Part B/Replot/Plots/To edit/qPCR")

PSAH12<-ggplot(data=FoldChange_stat %>% filter(target =="PsaH", timepoint ==12), 
               aes(x=strain, y=mean_lg2_FC,ymin= mean_lg2_FC-sd, ymax = mean_lg2_FC+sd))+
  geom_bar(stat="identity",position = "dodge",width=0.3, fill="grey35",color="black")+
  geom_errorbar(stat="identity",width = 0.15, color = "black")+
  geom_jitter(data=FoldChange %>%filter(target =="PsaH", timepoint ==12), 
              aes(x=strain,y=lg2FC),inherit.aes = FALSE, width = 0.08, size = 3, shape = 21, fill = "white")+
  scale_x_discrete(labels = xlabels)+
  scale_y_continuous(limits = c(-8,10), expand=c(0,0),breaks = seq(-8,10,2))+
  geom_hline(yintercept = 0,size=0.5) +
  ylab(bquote('Log'[2]~ 'Fold Change at 12 h'))+
  xlab(NULL)+
  theme_classic()+
  ggtitle(expression(paste(italic("PSAH"))))+
  theme(plot.title = element_text(hjust = 0.5,size =25))+
  theme(axis.title.y = element_text(size = 20),axis.text=element_text(size=16),axis.line.x = element_blank(),axis.ticks.x=element_blank())+
  theme(axis.text.x = element_text(angle = 40, hjust = 1))

PSAH12
ggsave(filename = "PSAH_qPCR_12.pdf", plot=last_plot(),width=7,height=8.5,units = "cm",scale=2,
       path="/Users/Meeko/Sync/Google Drive/Niyogi Lab/Zofingiensis/Switch Paper/Switch Paper Part B/Replot/Plots/To edit/qPCR")

PsbO12<-ggplot(data=FoldChange_stat %>% filter(target =="PsbO", timepoint ==12), 
               aes(x=strain, y=mean_lg2_FC,ymin= mean_lg2_FC-sd, ymax = mean_lg2_FC+sd))+
  geom_bar(stat="identity",position = "dodge",width=0.3, fill="grey35",color="black")+
  geom_errorbar(stat="identity",width = 0.15, color = "black")+
  geom_jitter(data=FoldChange %>%filter(target =="PsbO", timepoint ==12), 
              aes(x=strain,y=lg2FC),inherit.aes = FALSE, width = 0.08, size = 3, shape = 21, fill = "white")+
  scale_x_discrete(labels = xlabels)+
  scale_y_continuous(limits = c(-8,10), expand=c(0,0),breaks = seq(-8,10,2))+
  geom_hline(yintercept = 0,size=0.5) +
  ylab(bquote('Log'[2]~ 'Fold Change at 12 h'))+
  xlab(NULL)+
  theme_classic()+
  ggtitle(expression(paste(italic("PsbO"))))+
  theme(plot.title = element_text(hjust = 0.5,size =25))+
  theme(axis.title.y = element_text(size = 20),axis.text=element_text(size=16),axis.line.x = element_blank(),axis.ticks.x=element_blank())+
  theme(axis.text.x = element_text(angle = 40, hjust = 1))

PsbO12
ggsave(filename = "PsbO_qPCR_12.pdf", plot=last_plot(),width=7,height=8.5,units = "cm",scale=2,
       path="/Users/Meeko/Sync/Google Drive/Niyogi Lab/Zofingiensis/Switch Paper/Switch Paper Part B/Replot/Plots/To edit/qPCR")

RBCS112<-ggplot(data=FoldChange_stat %>% filter(target =="RBCS1", timepoint ==12), 
                aes(x=strain, y=mean_lg2_FC,ymin= mean_lg2_FC-sd, ymax = mean_lg2_FC+sd))+
  geom_bar(stat="identity",position = "dodge",width=0.3, fill="grey35",color="black")+
  geom_errorbar(stat="identity",width = 0.15, color = "black")+
  geom_jitter(data=FoldChange %>%filter(target =="RBCS1", timepoint ==12), 
              aes(x=strain,y=lg2FC),inherit.aes = FALSE, width = 0.08, size = 3, shape = 21, fill = "white")+
  scale_x_discrete(labels = xlabels)+
  scale_y_continuous(limits = c(-8,10), expand=c(0,0),breaks = seq(-8,10,2))+
  geom_hline(yintercept = 0,size=0.5) +
  ylab(bquote('Log'[2]~ 'Fold Change at 12 h'))+
  xlab(NULL)+
  theme_classic()+
  ggtitle(expression(paste(italic("RBCS1"))))+
  theme(plot.title = element_text(hjust = 0.5,size =25))+
  theme(axis.title.y = element_text(size = 20),axis.text=element_text(size=16),axis.line.x = element_blank(),axis.ticks.x=element_blank())+
  theme(axis.text.x = element_text(angle = 40, hjust = 1))

RBCS112
ggsave(filename = "RBCS1_qPCR_12.pdf", plot=last_plot(),width=7,height=8.5,units = "cm",scale=2,
       path="/Users/Meeko/Sync/Google Drive/Niyogi Lab/Zofingiensis/Switch Paper/Switch Paper Part B/Replot/Plots/To edit/qPCR")

MLDP112<-ggplot(data=FoldChange_stat %>% filter(target =="MLDP1", timepoint ==12), 
                aes(x=strain, y=mean_lg2_FC,ymin= mean_lg2_FC-sd, ymax = mean_lg2_FC+sd))+
  geom_bar(stat="identity",position = "dodge",width=0.3, fill="grey35",color="black")+
  geom_errorbar(stat="identity",width = 0.15, color = "black")+
  geom_jitter(data=FoldChange %>%filter(target =="MLDP1", timepoint ==12), 
              aes(x=strain,y=lg2FC),inherit.aes = FALSE, width = 0.08, size = 3, shape = 21, fill = "white")+
  scale_x_discrete(labels = xlabels)+
  scale_y_continuous(limits = c(-8,10), expand=c(0,0),breaks = seq(-8,10,2))+
  geom_hline(yintercept = 0,size=0.5) +
  ylab(bquote('Log'[2]~ 'Fold Change at 12 h'))+
  xlab(NULL)+
  theme_classic()+
  ggtitle(expression(paste(italic("MLDP1"))))+
  theme(plot.title = element_text(hjust = 0.5,size =25))+
  theme(axis.title.y = element_text(size = 20),axis.text=element_text(size=16),axis.line.x = element_blank(),axis.ticks.x=element_blank())+
  theme(axis.text.x = element_text(angle = 40, hjust = 1))

MLDP112
ggsave(filename = "MLDP1_qPCR_12.pdf", plot=last_plot(),width=7,height=8.5,units = "cm",scale=2,
       path="/Users/Meeko/Sync/Google Drive/Niyogi Lab/Zofingiensis/Switch Paper/Switch Paper Part B/Replot/Plots/To edit/qPCR")

SMP112<-ggplot(data=FoldChange_stat %>% filter(target =="SMP1", timepoint ==12), 
               aes(x=strain, y=mean_lg2_FC,ymin= mean_lg2_FC-sd, ymax = mean_lg2_FC+sd))+
  geom_bar(stat="identity",position = "dodge",width=0.3, fill="grey35",color="black")+
  geom_errorbar(stat="identity",width = 0.15, color = "black")+
  geom_jitter(data=FoldChange %>%filter(target =="SMP1", timepoint ==12), 
              aes(x=strain,y=lg2FC),inherit.aes = FALSE, width = 0.08, size = 3, shape = 21, fill = "white")+
  scale_x_discrete(labels = xlabels)+
  scale_y_continuous(limits = c(-8,10), expand=c(0,0),breaks = seq(-8,10,2))+
  geom_hline(yintercept = 0,size=0.5) +
  ylab(bquote('Log'[2]~ 'Fold Change at 12 h'))+
  xlab(NULL)+
  theme_classic()+
  ggtitle(expression(paste(italic("SMP1"))))+
  theme(plot.title = element_text(hjust = 0.5,size =25))+
  theme(axis.title.y = element_text(size = 20),axis.text=element_text(size=16),axis.line.x = element_blank(),axis.ticks.x=element_blank())+
  theme(axis.text.x = element_text(angle = 40, hjust = 1))

SMP112
ggsave(filename = "SMP1_qPCR_12.pdf", plot=last_plot(),width=7,height=8.5,units = "cm",scale=2,
       path="/Users/Meeko/Sync/Google Drive/Niyogi Lab/Zofingiensis/Switch Paper/Switch Paper Part B/Replot/Plots/To edit/qPCR")

PYK212<-ggplot(data=FoldChange_stat %>% filter(target =="PYK2", timepoint ==12), 
               aes(x=strain, y=mean_lg2_FC,ymin= mean_lg2_FC-sd, ymax = mean_lg2_FC+sd))+
  geom_bar(stat="identity",position = "dodge",width=0.3, fill="grey35",color="black")+
  geom_errorbar(stat="identity",width = 0.15, color = "black")+
  geom_jitter(data=FoldChange %>%filter(target =="PYK2", timepoint ==12), 
              aes(x=strain,y=lg2FC),inherit.aes = FALSE, width = 0.08, size = 3, shape = 21, fill = "white")+
  scale_x_discrete(labels = xlabels)+
  scale_y_continuous(limits = c(-8,10), expand=c(0,0),breaks = seq(-8,10,2))+
  geom_hline(yintercept = 0,size=0.5) +
  ylab(bquote('Log'[2]~ 'Fold Change at 12 h'))+
  xlab(NULL)+
  theme_classic()+
  ggtitle(expression(paste(italic("PYK2"))))+
  theme(plot.title = element_text(hjust = 0.5,size =25))+
  theme(axis.title.y = element_text(size = 20),axis.text=element_text(size=16),axis.line.x = element_blank(),axis.ticks.x=element_blank())+
  theme(axis.text.x = element_text(angle = 40, hjust = 1))

PYK212
ggsave(filename = "PYK2_qPCR_12.pdf", plot=last_plot(),width=7,height=8.5,units = "cm",scale=2,
       path="/Users/Meeko/Sync/Google Drive/Niyogi Lab/Zofingiensis/Switch Paper/Switch Paper Part B/Replot/Plots/To edit/qPCR")

LHC1612<-ggplot(data=FoldChange_stat %>% filter(target =="LHC16", timepoint ==12), 
                aes(x=strain, y=mean_lg2_FC,ymin= mean_lg2_FC-sd, ymax = mean_lg2_FC+sd))+
  geom_bar(stat="identity",position = "dodge",width=0.3, fill="grey35",color="black")+
  geom_errorbar(stat="identity",width = 0.15, color = "black")+
  geom_jitter(data=FoldChange %>%filter(target =="LHC16", timepoint ==12), 
              aes(x=strain,y=lg2FC),inherit.aes = FALSE, width = 0.08, size = 3, shape = 21, fill = "white")+
  scale_x_discrete(labels = xlabels)+
  scale_y_continuous(limits = c(-8,10), expand=c(0,0),breaks = seq(-8,10,2))+
  geom_hline(yintercept = 0,size=0.5) +
  ylab(bquote('Log'[2]~ 'Fold Change at 12 h'))+
  xlab(NULL)+
  theme_classic()+
  ggtitle(expression(paste(italic("LHC16"))))+
  theme(plot.title = element_text(hjust = 0.5,size =25))+
  theme(axis.title.y = element_text(size = 20),axis.text=element_text(size=16),axis.line.x = element_blank(),axis.ticks.x=element_blank())+
  theme(axis.text.x = element_text(angle = 40, hjust = 1))

LHC1612
ggsave(filename = "LHC16_qPCR_12.pdf", plot=last_plot(),width=7,height=8.5,units = "cm",scale=2,
       path="/Users/Meeko/Sync/Google Drive/Niyogi Lab/Zofingiensis/Switch Paper/Switch Paper Part B/Replot/Plots/To edit/qPCR")



COX1012<-ggplot(data=FoldChange_stat %>% filter(target =="COX10", timepoint ==12), 
                aes(x=strain, y=mean_lg2_FC,ymin= mean_lg2_FC-sd, ymax = mean_lg2_FC+sd))+
  geom_bar(stat="identity",position = "dodge",width=0.3, fill="grey35",color="black")+
  geom_errorbar(stat="identity",width = 0.15, color = "black")+
  geom_jitter(data=FoldChange %>%filter(target =="COX10", timepoint ==12), 
              aes(x=strain,y=lg2FC),inherit.aes = FALSE, width = 0.08, size = 3, shape = 21, fill = "white")+
  scale_x_discrete(labels = xlabels)+
  scale_y_continuous(limits = c(-8,10), expand=c(0,0),breaks = seq(-8,10,2))+
  geom_hline(yintercept = 0,size=0.5) +
  ylab(bquote('Log'[2]~ 'Fold Change at 12 h'))+
  xlab(NULL)+
  theme_classic()+
  ggtitle(expression(paste(italic("COX10"))))+
  theme(plot.title = element_text(hjust = 0.5,size =25))+
  theme(axis.title.y = element_text(size = 20),axis.text=element_text(size=16),axis.line.x = element_blank(),axis.ticks.x=element_blank())+
  theme(axis.text.x = element_text(angle = 40, hjust = 1))

COX1012
ggsave(filename = "COX10_qPCR_12.pdf", plot=last_plot(),width=7,height=8.5,units = "cm",scale=2,
       path="/Users/Meeko/Sync/Google Drive/Niyogi Lab/Zofingiensis/Switch Paper/Switch Paper Part B/Replot/Plots/To edit/qPCR")

SAD112<-ggplot(data=FoldChange_stat %>% filter(target =="SAD1", timepoint ==12), 
               aes(x=strain, y=mean_lg2_FC,ymin= mean_lg2_FC-sd, ymax = mean_lg2_FC+sd))+
  geom_bar(stat="identity",position = "dodge",width=0.3, fill="grey35",color="black")+
  geom_errorbar(stat="identity",width = 0.15, color = "black")+
  geom_jitter(data=FoldChange %>%filter(target =="SAD1", timepoint ==12), 
              aes(x=strain,y=lg2FC),inherit.aes = FALSE, width = 0.08, size = 3, shape = 21, fill = "white")+
  scale_x_discrete(labels = xlabels)+
  scale_y_continuous(limits = c(-8,10), expand=c(0,0),breaks = seq(-8,10,2))+
  geom_hline(yintercept = 0,size=0.5) +
  ylab(bquote('Log'[2]~ 'Fold Change at 12 h'))+
  xlab(NULL)+
  theme_classic()+
  ggtitle(expression(paste(italic("SAD1"))))+
  theme(plot.title = element_text(hjust = 0.5,size =25))+
  theme(axis.title.y = element_text(size = 20),axis.text=element_text(size=16),axis.line.x = element_blank(),axis.ticks.x=element_blank())+
  theme(axis.text.x = element_text(angle = 40, hjust = 1))

SAD112
ggsave(filename = "SAD1_qPCR_12.pdf", plot=last_plot(),width=7,height=8.5,units = "cm",scale=2,
       path="/Users/Meeko/Sync/Google Drive/Niyogi Lab/Zofingiensis/Switch Paper/Switch Paper Part B/Replot/Plots/To edit/qPCR")

dev.off()
