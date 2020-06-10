###############################################################
#AUTHOR:Shivani
#This is a script for plotting heatmaps using the RNA_Seq masterdataset and the logFC values calculated using the RNA Seq data

################################################################

#loading the useful libraries#

library(tidyverse)
library(ggplot2)
library(pheatmap)

#reading the transcription factor file containing the geneIDs

upregulated_txn_tf2 =  separate(readxl::read_xlsx("upregulated_txn.xlsx"), Description, c("geneID", "NA","add"), sep = " -")


#But, need to make a list with FC and GO terms()

RNA_logFC = readxl::read_xlsx("logFCvalues.xlsx")

#reading the masterdataset files with FPKM values
RNA_master = readxl::read_xlsx("RNA_masterdataset.xlsx")


#this is to extract logFC values for only upregulated_txn from the total logFC file

#this is to extract logFC values for only myb from the total logFC file

##########################################################
#YOU ONLY MAKE CHANGES IN THE FOLLOWING LINE TO GET THE REQUIRED FOLD CHANGE CUTOFF!
###########################################################
logFC_upregulated_txn = filter(RNA_logFC, geneID %in% (upregulated_txn_tf2 %>% pull(geneID))) %>% filter_if(is.numeric, any_vars(abs(.)>=2)) %>% replace(., is.na(.), 0)


########################################################################

fpkmFC_upregulated_txn = filter(RNA_master, geneID %in% (logFC_upregulated_txn %>% pull(geneID))) %>% replace(., is.na(.), 0) 


fpkmFC_upregulated_txn_cutoff = fpkmFC_upregulated_txn %>% filter_if(is.numeric, all_vars(abs(.)>=1)) %>% replace(., is.na(.), 0)

logFC_upregulated_txn_cutoff = filter(logFC_upregulated_txn, geneID %in% (fpkmFC_upregulated_txn_cutoff %>% pull(geneID)))

fpkm_log2 = readxl::read_xlsx("log2_fpkm.xlsx")


##################################################################################################################
#for plotting FPKM


fpkm_log2 = readxl::read_xlsx("log2_fpkm.xlsx")

log2fpkm_upregulated_txn = filter(fpkm_log2, geneID %in% (logFC_upregulated_txn_cutoff %>% pull(geneID))) %>% replace(., is.na(.), 0)

log2fpkm_upregulated_txn %>% column_to_rownames("geneID") -> heatmap_data_FPKM

logFC_upregulated_txn_cutoff %>% column_to_rownames("geneID") -> heatmap_data_logFC



##################################################################################################################
#annotations for the heatmap



anntn =  readxl::read_xlsx("annotations.xlsx")

##################################################################################################################
#predictions using DeepLoc
#although the list says predalgo, the algorithm used is deeploc
predalgo = separate(readxl::read_xlsx("predAlgo.xlsx"), ID, c("geneID", "NA","extra"), sep = ".t1")

#final annotations file

upregulated_txn_anntn = upregulated_txn_tf2 %>% left_join(anntn, by = "geneID") %>% left_join(predalgo, by = "geneID") %>% replace(., is.na(.), 0)


upregulated_txn_anntn %>% dplyr::select(geneID,GO_terms, Target) %>% column_to_rownames("geneID") -> row_annotations

#getting all the annotations together

data.frame(row.names = colnames(heatmap_data_logFC), 
           Treatment = factor(grepl("-G$",colnames(heatmap_data_logFC)),
                              levels = c(T,F),labels = c("Remove Glu","Add Glu"))) -> col_anntn


#scaling the heatmap

myScale4x4<-c(seq(-4,-1.5,length=60),
              seq(-1.499,-0.5,length=60),
              seq(-0.499,0.499,length=60),
              seq(0.5,1.499,length=60),
              seq(1.5,4, length=60))
myScale3x3<-c(seq(-3,-1.5,length=60),
              seq(-1.499,-0.5,length=60),
              seq(-0.499,0.499,length=60),
              seq(0.5,1.499,length=60),
              seq(1.5,3, length=60))

#using standardized colors for the heatmap
bbyColors<-colorRampPalette(c('dodgerblue','black','yellow'))(n=299)
#scaling the output file
pdfheight = (nrow(heatmap_data_logFC)*0.5+2)
pdfwidth = (ncol(heatmap_data_logFC)*0.5+10)


#finally plotting the heatmap
heatmap_data_logFC %>% pheatmap(cluster_cols = F, color=bbyColors,
                                cluster_rows = F,
                                breaks = myScale4x4,
                                #kmeans_k = 2,
                                annotation_row = row_annotations, 
                                annotation_col = col_anntn,
                                #annotation_colors = mycolors,
                                annotation_colors = c(list(Treatment=c("Add Glu" = "orange","Remove Glu" = "green")),GO_terms = bbyColors),
                                angle_col = 45, 
                                treeheight_row = 10,
                                show_rownames = T,
                                display_numbers = T,
                                gaps_col = 5,
                                filename = "upregulated_txn_heatmap.pdf", width = pdfwidth, height = pdfheight,
                                cellwidth = 40, cellheight = 15
                                
)

#defining colors for plotting FPKM heatmap
fmColors<-colorRampPalette(c("yellow","red","brown2"))(n=299)


myScale0x20<-c(seq(0,5,length=60),
               seq(4.99,8,length=60),
               seq(7.99,12,length=60),
               seq(11.99,16,length=60),
               seq(15.99,20, length=60))


#FPKM heatmap
heatmap_data_FPKM %>% pheatmap(cluster_cols = F, color=fmColors,
                               cluster_rows = F,
                               
                               breaks = myScale0x20,
                               angle_col = 45, 
                               #kmeans_k = 2,
                               show_rownames = T,
                               #cutree_cols = 4,
                               gaps_col = c(6,11,17),
                               display_numbers = T,
                               filename = "upregulated_txn_heatmap_fpkm.pdf", width = pdfwidth, height = pdfheight,
                               cellwidth = 40, cellheight = 15
                               
)







