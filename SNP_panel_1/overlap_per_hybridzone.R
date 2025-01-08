### compare 99.9perc Fsts from 3 different datasets: fstc, fstv & angsd for autosomes and chrz separately ###
### Rscript overlap_per_hybridzone.R G01_G02 G06_G07 G08_G09 chr_autosomes chr_chrz ###

#Rversion 4.0.2

library(dplyr)
library(devtools)
library(grid)
library(ggplot2)
library(ggvenn)
library(VennDiagram)
library(base)

args <- commandArgs(trailingOnly = TRUE)

for (chrtype in args[4:5]) 
{
cat(chrtype) 
for (filename in args[1:3]) 
{
cat(filename)
setwd("/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/01_probes/fst/fst_summary")
datc <- read.delim(file=paste("fst_christen",filename,chrtype,"only.txt",sep="_"), header=TRUE, sep="\t")
fstc <- datc$WEIR_AND_COCKERHAM_FST
nrow(na.omit(datc))
fstc[fstc<0] <- 0 #converting all negative fsts to 0
mean(fstc, na.rm=TRUE) 
quantile(fstc, c(.80, .95, .99, .999, .9999, .99995, .99999), na.rm=TRUE)
datc99.9 <- datc %>% filter(fstc > quantile(fstc, 0.999, na.rm=TRUE))
nrow(datc99.9)

datv <- read.delim(file=paste("fst_verena",filename,chrtype,"only.txt",sep="_"), header=TRUE, sep="\t")
fstv <- datv$WEIR_AND_COCKERHAM_FST
nrow(na.omit(datv)) 
fstv[fstv<0] <- 0
mean(fstv, na.rm=TRUE)
quantile(fstc, c(.80, .95, .99, .999, .9999, .99995, .99999), na.rm=TRUE)
datv99.9 <- datv %>% filter(fstv > quantile(fstv, 0.999, na.rm=TRUE))
nrow(datv99.9)

datangsd <- read.delim(file=paste("fst_angsd",filename,chrtype,"only.txt",sep="_"), header=TRUE, sep="\t")
fstangsd <- datangsd$FST
nrow(na.omit(datangsd)) #no NA or neg values in angsd
mean(fstangsd, na.rm=TRUE) 
quantile(fstangsd, c(.80, .95, .99, .999, .9999, .99995, .99999), na.rm=TRUE)
datangsd99.9 <- datangsd %>% filter(fstangsd > quantile(fstangsd, 0.999, na.rm=TRUE))
nrow(datangsd99.9) 

#check overlap
xlist <- list(GATK_Christen=datc99.9$LocusName, GATK_Verena=datv99.9$LocusName, ANGSD=datangsd99.9$LocusName)
x <- c(datc99.9$LocusName, datv99.9$LocusName, datangsd99.9$LocusName)
y <- data.frame(unique(x, incomparables = FALSE))
nrow(y) 
pdf(file=paste(filename,chrtype,"fst99.9",sep="_"))
ggvenn(xlist, 
       fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
       stroke_size = 0.5, set_name_size = 4) +
       labs(title=paste(filename,chrtype,"99.9",sep="_")) + 
       theme(plot.title = element_text(hjust=0.5))
dev.off()
overlap <- calculate.overlap(xlist)
summary(overlap)
overlap_c <- data.frame(LocusName=overlap$a1, OverlapType="C")
overlap_cv <- data.frame(LocusName=overlap$a2, OverlapType="CV")
overlap_v <- data.frame(LocusName=overlap$a3, OverlapType="V")
overlap_ca <- data.frame(LocusName=overlap$a4, OverlapType="CA")
overlap_cva <- data.frame(LocusName=overlap$a5, OverlapType="CVA")
overlap_va <- data.frame(LocusName=overlap$a6, OverlapType="VA")
overlap_a <- data.frame(LocusName=overlap$a7, OverlapType="A")
overlap_summary <- merge(overlap_a, merge(overlap_c, merge(overlap_v, merge(overlap_va, merge(overlap_ca, merge(overlap_cva, overlap_cv, all=TRUE), all=TRUE), all=TRUE), all=TRUE), all=TRUE), all=TRUE)
overlap_summary$FstChristen <- overlap_summary$LocusName #dummy
overlap_summary$FstVerena <- overlap_summary$LocusName #dummy
overlap_summary$FstAngsd <- overlap_summary$LocusName #dummy
write.table(overlap_summary, file=paste(filename,chrtype,"persite_fst99.9.txt",sep="_"), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
}
}

