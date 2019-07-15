##don't forget the library dependencies
require(dplyr)
require(ggplot2)

#For these files, I used a separate text file to manually seek and replace gene names and time points in order to generate all 28 or so plots. Then I faceted them in Illustrator. 

png("/home/daniel/Sync/Google Drive/Niyogi Lab/Zofingiensis/Switch Paper/qPCR/qPCR rPlots/gene.png", height = 600, width = 600)
qpcr <- read.csv("/home/daniel/Sync/Google Drive/Niyogi Lab/Zofingiensis/Switch Paper/qPCR/qPCR_Data_all2.csv", header=TRUE,as.is=TRUE) %>% tbl_df
qpcr <- qpcr %>% group_by(strain, bio_reps, treatment, timepoint, rxn, target) %>% mutate(mean_ct = mean(ct)) %>% ungroup
qpcr_apps <- qpcr %>% filter(target == "APPS") %>% group_by(strain, bio_reps, tech_reps, treatment, timepoint, rxn) %>% summarize(apps_ct = mean(mean_ct))
qpcr2 <- left_join(qpcr, qpcr_apps, by=c("strain", "bio_reps", "tech_reps", "treatment", "timepoint", "rxn")) %>% filter(target != "APPS")
qpcr3 <- qpcr2 %>% mutate(delta_ct = mean_ct - apps_ct)
tmp <- qpcr3 %>% distinct(strain, treatment, target, timepoint, bio_reps, delta_ct)


p1 <- tmp %>% filter(target == "HXK1", timepoint == 0.5) %>% ggplot(aes( x= interaction(treatment,strain),y=delta_ct)) + 
  geom_jitter(width=0.001, alpha=1, color = "blue") +
  geom_boxplot(alpha = 0.1)+
  ggtitle(expression(paste(italic("HXK1")," 0.5h"))) + 
  xlab("") +
  ylab("Delta CT") +
  theme(plot.title = element_text(color = "black", size = 14, face = "bold", hjust = 0.5))+theme_bw(base_size = 20) + ylim(-13,13)
p1
dev.off()
#ggsave(p1, filename="HXK_12h")


