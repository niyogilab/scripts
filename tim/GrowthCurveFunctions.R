#Functions for Zofingiensis growth curves

#Make sure you setwd to where you store your csv files. Go to Session > "Set Working Directory" > "Choose Directory"
#and pick the file that your csv files are saved

#Make sure you have ggplot2 and aggregate installed
install.packages("ggplot2")
library(ggplot2)

#activate all ofthe following functions. R will save the function after you activatein a session. 
standard.error <- function(x){
  sd(x)/sqrt(length(x))
}

densitybyday<-function(x){
  file<-read.csv(x, header=T)
  agg<-do.call(data.frame, 
               aggregate(cbind(file$CellDensity,file$Mean,file$VOL) ~ file$Media + file$Day,
                         data=file,
                         FUN=function(file) c(mn=mean(file), se=standard.error(file))))
  colnames(agg) <- c("media", "day", "celldensity","celldensity.se", "size", "size.se", "totalvol", "totalvol.se")
  celldensity <- ggplot(agg, aes(x=day, y=celldensity, colour=media))
  celldensity + geom_line() + geom_point()+geom_errorbar(aes(ymin=celldensity-celldensity.se, ymax=celldensity+celldensity.se),width=0.1) +  xlim(0, max(agg$day))
   }

volbyday<-function(x){
  file<-read.csv(x, header=T)
  agg<-do.call(data.frame, 
               aggregate(cbind(file$CellDensity,file$Mean,file$TotalVol) ~ file$Media + file$Day,
                         data=file,
                         FUN=function(file) c(mn=mean(file), se=standard.error(file))))
  colnames(agg) <- c("media", "day", "celldensity","celldensity.se", "size", "size.se", "totalvol", "totalvol.se")
  celldensity <- ggplot(agg, aes(x=day, y=totalvol, colour=media))
  celldensity + geom_line() + geom_point()+geom_errorbar(aes(ymin=totalvol-totalvol.se, ymax=totalvol+totalvol.se),width=0.1) + xlim(0, max(agg$day))
}

meansizebyday<-function(x){
  file<-read.csv(x, header=T)
  agg<-do.call(data.frame, 
               aggregate(cbind(file$CellDensity,file$Mean,file$VOL) ~ file$Media + file$Day,
                         data=file,
                         FUN=function(file) c(mn=mean(file), se=standard.error(file))))
  colnames(agg) <- c("media", "day", "celldensity","celldensity.se", "meansize", "meansize.se", "totalvol", "totalvol.se")
  celldensity <- ggplot(agg, aes(x=day, y=meansize, colour=media))
  celldensity + geom_line() + geom_point()+geom_errorbar(aes(ymin=meansize-meansize.se, ymax=meansize+meansize.se),width=0.1) + xlim(0, max(agg$day)) + ylim(2, 10)
}

P1byday<-function(x){
  file<-read.csv(x, header=T)
  agg<-do.call(data.frame, 
               aggregate(cbind(file$CellDensity,file$P1,file$VOL) ~ file$Media + file$Day,
                         data=file,
                         FUN=function(file) c(mn=mean(file), se=standard.error(file))))
  colnames(agg) <- c("media", "day", "celldensity","celldensity.se", "P1", "P1.se", "totalvol", "totalvol.se")
  celldensity <- ggplot(agg, aes(x=day, y=P1, colour=media))
  celldensity + geom_line() + geom_point()+geom_errorbar(aes(ymin=P1-P1.se, ymax=P1+P1.se),width=0.1)  + xlim(0, max(agg$day)) + ylim(2, 10)
}

FvFmbyday<-function(x){
  file<-read.csv(x, header=T)
  agg<-do.call(data.frame, 
               aggregate(cbind(file$CellDensity,file$Fv.Fm,file$VOL) ~ file$Media + file$Day,
                         data=file,
                         FUN=function(file) c(mn=mean(file), se=standard.error(file))))
  colnames(agg) <- c("media", "day", "celldensity","celldensity.se", "Fv.Fm", "Fv.Fm.se", "totalvol", "totalvol.se")
  celldensity <- ggplot(agg, aes(x=day, y=Fv.Fm, colour=media))
  celldensity + geom_line() + geom_point()+geom_errorbar(aes(ymin=Fv.Fm-Fv.Fm.se, ymax=Fv.Fm+Fv.Fm.se),width=0.1)  + xlim(0, max(agg$day)) + ylim(0, 1)
}


pHbyday<-function(x){
  file<-read.csv(x, header=T)
  agg<-do.call(data.frame, 
               aggregate(cbind(file$CellDensity,file$pH,file$VOL) ~ file$Media + file$Day,
                         data=file,
                         FUN=function(file) c(mn=mean(file), se=standard.error(file))))
  colnames(agg) <- c("media", "day", "celldensity","celldensity.se", "pH", "pH.se", "totalvol", "totalvol.se")
  pH <- ggplot(agg, aes(x=day, y=pH, colour=media))
  pH + geom_line() + geom_point()+geom_errorbar(aes(ymin=pH-pH.se, ymax=pH+pH.se),width=0.1) + xlim(0, max(agg$day))
}

peaknumber <- function(x){
  file<-read.csv(x, header=T)
  agg<-do.call(data.frame, 
               aggregate(cbind(file$CellDensity,file$PeakNum,file$VOL) ~ file$Media + file$Day,
                         data=file,
                         FUN=function(file) c(mn=mean(file), se=standard.error(file))))
  colnames(agg) <- c("media", "day", "celldensity","celldensity.se", "PeakNum", "PeakNum.se", "totalvol", "totalvol.se")
  PeakNum <- ggplot(agg, aes(x=day, y=PeakNum, colour=media))
  PeakNum + geom_line() + geom_point()+geom_errorbar(aes(ymin=PeakNum-PeakNum.se, ymax=PeakNum+PeakNum.se),width=0.1) + xlim(0, max(agg$day))
}



#Pipeline for seeing all of the parameters
file <- "yourfile.csv"


#functions that will map trait each day. Input your file: 

densitybyday(file) + ylab(expression(paste("Cells ", mL^{-1})))
volbyday(file) + ylab(expression(paste(mu*m^{3}, " ", mL^{-1})))
meansizebyday(file) 
P1byday(file)
FvFmbyday(file) 
pHbyday(file)  


#copy and paste graphs by going to export -> "copy to clipboard" then pasting onto a word doc, etc. 



file <-read.csv("FvFmOff.csv", header=T)
FvFm <- do.call(data.frame, aggregate(file$FvFm ~ file$Sample + file$Time,
          data=file,
          FUN=function(file) c(mn=mean(file), se=standard.error(file))))
colnames(FvFm) <- c("Sample", "Time", "FvFm", "FvFm.se")
o <- ggplot(FvFm, aes(x=Time, y=FvFm, colour=Sample, ymin=FvFm-FvFm.se, ymax=FvFm+FvFm.se)) 
o + geom_line(size=0.8) +geom_errorbar() + geom_point()+ geom_errorbar(width=0.1) + theme_bw()  

#Let's do it just 5 hours
#FvFm <- FvFm[which(FvFm$Time<15),]
o <- ggplot(FvFm, aes(x=Time, y=FvFm, colour=Sample, ymin=FvFm-FvFm.se, ymax=FvFm+FvFm.se)) 
o + geom_line(size=0.8) +geom_errorbar() + geom_point()+ geom_errorbar(width=0.1) + theme_bw() + xlim(0,5)



file <-read.csv("FvFmOff2.csv", header=T)
FvFm <- do.call(data.frame, aggregate(file$FvFm ~ file$Sample + file$Time,
                                      data=file,
                                      FUN=function(file) c(mn=mean(file), se=standard.error(file))))
colnames(FvFm) <- c("Sample", "Time", "FvFm", "FvFm.se")
o <- ggplot(FvFm, aes(x=Time, y=FvFm, colour=Sample, ymin=FvFm-FvFm.se, ymax=FvFm+FvFm.se)) 
o + geom_line(size=0.8) +geom_errorbar() + geom_point()+ geom_errorbar(width=0.1) + theme_bw() + 
  
#Just get the wt and glucose
c(which(FvFm$Sample[1]==factor(Wt)),  which(FvFm$Sample== "Wt +Glc+DMSO"))

FvFm <- FvFm[c(which(as.character(FvFm$Sample)=="Wt "),which(as.character(FvFm$Sample)=="Wt +Glc+DMSO")),] 
o <- ggplot(FvFm, aes(x=Time, y=FvFm, colour=Sample, ymin=FvFm-FvFm.se, ymax=FvFm+FvFm.se)) 
o + geom_line(size=0.8) +geom_errorbar() + geom_point()+ geom_errorbar(width=0.1) + theme_bw()

file <-read.csv("FvFmOff2.csv", header=T)
FvFm <- do.call(data.frame, aggregate(file$FvFm ~ file$Sample + file$Time,
                                      data=file,
                                      FUN=function(file) c(mn=mean(file), se=standard.error(file))))
colnames(FvFm) <- c("Sample", "Time", "FvFm", "FvFm.se")
FvFm <- FvFm[-which(as.character(FvFm$Sample)=="Wt +Glc+CHX"),]
o <- ggplot(FvFm, aes(x=Time, y=FvFm, colour=Sample, ymin=FvFm-FvFm.se, ymax=FvFm+FvFm.se)) 
o + geom_line(size=0.8) +geom_errorbar() + geom_point()+ geom_errorbar(width=0.1) + theme_bw()


