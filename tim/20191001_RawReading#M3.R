#Reading M3 Files
#This code simply extracts BinDiam and BinHeight from raw files that you can then plot with the filter function
#You want to keep all your #M3 files in a single folder that you save in the working directory here this is happening

#inputs: 

# The working directory that contains the folder with files .#M3 is here should be in path:
path = ""
#input the folder 
folder= ""
  
#run the code below and for loop. If you use average files: sorry they are evil and I don't know how to process them.
  #It's not good to use the average function on the coulter counter anyway. just do technical replicates and average later. 

#this will make a list of the files 
filelist <- list.files(x) ;filelist<- filelist[!grepl("av|pm|starte", filelist)] ;
#Set up the counterdata vector to get the files 
counterdata <- NULL
for(i in 1:length(filelist)){
#here is the rawest example of how to get the data
filename = filelist[i]
counterdata <- rbind(counterdata, cell.distrib)
Practice <- read.csv(file=paste(x,"/",filelist[i], sep=""),
                     row.names=NULL) #row.names=NULL because we don't have row names 
Practice <- Practice %>% mutate(RowNum = as.numeric(row.names(Practice))) #I want to get a numeric list of rows so I can find in betweens 

#The following steps extract the data related to diameter of cells BinDIAM and the number in those (BinHeight)
Points <- filter(Practice, row.names=="[#Bindiam]"| row.names=="[#Binheight]"|row.names=="[Binunits]" |row.names=="[SizeStats]")$RowNum  #This lists the values where the numeric values start and stop.
#I have a feeling they will alway be 490, 748, 751, 1009 
BinDiam = as.numeric(Practice$row.names[(Points[1]+1):(Points[2]-1)])
BinHeight = as.numeric(Practice$row.names[(Points[3]+1):(Points[4]-1)])  #The first 3 are like this 
#Potential Debugging point: are the lengths the same: 
length(BinDiam); length(BinHeight)
#Let's make a file with these two things together 
cell.distrib <- data_frame(filename, BinDiam, BinHeight) %>% filter(BinDiam>2) #don't need to add dilution because it is already multiplied in there
counterdata <- rbind(counterdata, cell.distrib)
}

#all size distributions are now in folder counterdata. It's a long file

#example plot, using geom_smooth. filelist specifies any filename you want to look at. you sould mess with the IDs and do faceting as well
# as time ordered for loops to get a better sense of organization 
counterdata %>% filter(filename == filelist[20]) %>% filter(BinDiam <20)  %>%  #you will want to mess with the BinDiam 
  ggplot() +
  stat_smooth(method="loess", span=0.1, se=TRUE, aes_string(x="BinDiam", y="BinHeight", group="filename", color="filename", fill="filename"), alpha=0.3) + theme_classic(base_size = 22) +
   ylab("Count") + xlab(expression(paste("Cell Diameter ", mu*m))) 




