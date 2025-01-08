### consolidate Fst outliers of all three hybridzones ###
### Rscript overlap_all_hybridzones.R G01_G02 G06_G07 G08_G09 ###

library(dplyr)
library(ggpubr)
library(devtools)
library(grid)
library(ggplot2)
library(ggvenn)
library(VennDiagram)
library(base)

setwd("/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/01_probes/fst/fst_summary")
args <- commandArgs(trailingOnly = TRUE)

#################### Further autosomal Fst filter ####################
for (filename in args[1:3]) 
{
cat(filename)
datc <- read.delim(file=paste("fst_christen",filename,"chr_autosomes_only.txt",sep="_"), header=TRUE, sep="\t")
datv <- read.delim(file=paste("fst_verena",filename,"chr_autosomes_only.txt",sep="_"), header=TRUE, sep="\t")
datangsd <- read.delim(file=paste("fst_angsd",filename,"chr_autosomes_only.txt",sep="_"), header=TRUE, sep="\t")
fstc[fstc<0] <- 0 
fstv[fstv<0] <- 0 
fstangsd[fstangsd<0] <- 0 

dat1 <- read.delim(file=paste(filename,"chr_autosomes_persite_fst99.9.txt",sep="_"), header=TRUE, sep="")
dat1[is.na(dat1)] <- 0
dat1[dat1<0] <- 0
dat1$t1 <- ifelse(dat1$FstChristen > quantile(fstc, 0.999,na.rm=TRUE) & dat1$FstVerena > 0 & dat1$FstAngsd > 0, "yes", "no")
dat1$t2 <- ifelse(dat1$FstVerena >  quantile(fstv, 0.999,na.rm=TRUE) & dat1$FstChristen > 0 & dat1$FstAngsd > 0, "yes", "no")
dat1$t3 <-ifelse(dat1$FstAngsd >  quantile(fstangsd, 0.999,na.rm=TRUE) & dat1$FstVerena > 0 & dat1$FstChristen > 0, "yes", "no")
dat1$t4 <-ifelse(dat1$FstAngsd > 0.9 | dat1$FstVerena > 0.9 | dat1$FstChristen > 0.9, "yes", "no")
dat1filtered <- dat1[apply(dat1, 1, function(v) any(grepl("yes", v))),]
nrow(dat1filtered)
write.table(dat1filtered, file=paste(filename,"persite_fst99.9_filtered.txt",sep="_"), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
}

#################### Outlier autosomal Fsts for all hybridz zones ####################
dat1filtered <- read.delim(file=paste(args[1],"persite_fst99.9_filtered.txt",sep="_"), sep="\t", row.names=FALSE, col.names=TRUE)
dat1filtered <- read.delim(file=paste(args[2],"persite_fst99.9_filtered.txt",sep="_"), sep="\t", row.names=FALSE, col.names=TRUE)
dat1filtered <- read.delim(file=paste(args[3],"persite_fst99.9_filtered.txt",sep="_"), sep="\t", row.names=FALSE, col.names=TRUE)
xlist <- list(G01vsG02=dat1filtered$LocusName, G06vsG07=dat2filtered$LocusName, G08vsG09= dat3filtered$LocusName)
overlap <- calculate.overlap(xlist)
summary(overlap)
overlap_G12 <- data.frame(LocusName=overlap$a1, OverlapType="G12")
overlap_G12_67 <- data.frame(LocusName=overlap$a2, OverlapType="G12_67")
overlap_G67 <- data.frame(LocusName=overlap$a3, OverlapType="G67")
overlap_G12_89 <- data.frame(LocusName=overlap$a4, OverlapType="G12_89")
overlap_G67_89 <- data.frame(LocusName=overlap$a6, OverlapType="G67_89")
overlap_G89 <- data.frame(LocusName=overlap$a7, OverlapType="G89")
overlap_summary <- merge(overlap_G89, merge(overlap_G67, merge(overlap_G12, merge(overlap_G67_89, merge(overlap_G12_67, overlap_G12_89, all=TRUE), all=TRUE), all=TRUE), all=TRUE), all=TRUE)
allhz <- as.data.frame(do.call(rbind, strsplit(allhz$LocusName, split=":")))
allhz$Scaffold <- overlap_summary_all$V1
allhz$Pos <- overlap_summary_all$V2
write.table(allhz, file="allhybridzones_persite_fst99.9_filtered.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
nrow(allhz)

#################### Outlier Fsts for all hybrid zones including chrz ####################
chrz <- read.delim("allhybridzones_chrz_fst99.9.txt", header=TRUE, sep="")
auto_chrz <- merge(allhz, chrz, all=TRUE)                                             
write.table(auto_chrz, file="./snp_panel_summary/allhybridzones_persite_chrzinc_fst99.9_filtered.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
nrow(auto_chrz)

#################### Outlier dxy for G01vsG02 ####################
setwd("/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/01_probes/dxy")
datG12 <- read.delim("dxypersite_G01.G02_christen.csv.nl.chr_autosomes_only.txt", header=TRUE, sep="\t")
datG12[is.na(datG12)] <- 0
datG12[datG12<0] <- 0
dxyG12 <- datG12$dxy_G01_G02
quantile(dxyG12, c(.80, .95, .99, .999, .9995, .9999, .99995, .99999), na.rm=TRUE) 
datG12_99.95 <- datG12 %>% filter(dxyG12 > quantile(dxyG12, 0.9995, na.rm=TRUE))
datG12_99.95$LocusName <- paste(datG12_99.95$scaffold,datG12_99.95$mid,sep=":")
nrow(datG12_99.95)       
q <- datG12_99.95[c(2,3,9,10,11,12)]
colnames(q) <- c("CHROM", "Pos", "Dxy", "Fst", "Scaffold", "LocusName")
write.table(q, file="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/01_probes/snp_panel_summary/G01_G02_persite_dxy99.95.txt"), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

#################### all_chrz_fst_dxy99.95_knief_strict ####################
setwd("/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/01_probes/snp_panel_summary)

kniefprobes <- read.delim("/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/01_probes/goldengate_snps/data_TableS1.txt", header=TRUE, sep="\t")
fstkniefchrzdxy <- list(fst99.9=as.character(allhz$LocusName), chrz=as.character(chrz$LocusName), dxy99.95=as.character(dxy99.95$LocusName), Kniefprobes=as.character(kniefprobes$LocusName))

pdf("all_chrz_fst_dxy99.95_knief_strict.pdf")
ggvenn(fstkniefchrzdxy, 
fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
stroke_size = 0.5, set_name_size = 3) +
labs(title="all_outier_sitesfiltereddxy99.95") + theme(plot.title = element_text(hjust=0.5))
dev.off()

overlap <- calculate.overlap(fstkniefchrzdxy)
summary(overlap)
overlap_a5 <- data.frame(LocusName=overlap$a5, OverlapType="fst_dxy_knief") 
overlap_a4 <- data.frame(LocusName=overlap$a4, OverlapType="fst_dxy")
overlap_a10 <- data.frame(LocusName=overlap$a10, OverlapType="fst_knief")
overlap_a8 <- data.frame(LocusName=overlap$a8, OverlapType="chrz_knief")
overlap_a2 <- data.frame(LocusName=overlap$a2, OverlapType="dxy_knief")
overlap_a9 <- data.frame(LocusName=overlap$a9, OverlapType="fst")
overlap_a14 <- data.frame(LocusName=overlap$a14, OverlapType="chrz")
overlap_a1 <- data.frame(LocusName=overlap$a1, OverlapType="dxy")
overlap_a3 <- data.frame(LocusName=overlap$a3, OverlapType="knief")
overlap_summary <- merge(overlap_a8, merge(overlap_a10, merge(overlap_a4, merge(overlap_a5, merge(overlap_a2,merge(overlap_a9, merge(overlap_a14, merge(overlap_a1, overlap_a3, all=TRUE), all=TRUE), all=TRUE), all=TRUE), all=TRUE), all=TRUE), all=TRUE), all=TRUE)
overlap_summary_all <- as.data.frame(do.call(rbind, strsplit(overlap_summary$LocusName, split=":")))
overlap_summary$Scaffold <- overlap_summary_all$V1
overlap_summary$Pos <- overlap_summary_all$V2
write.table(overlap_summary, file="./01_outlier/all_chrz_fst_dxy99.95_knief_strict_31k.txt", sep="\t", row.names=FALSE, col.names=TRUE)
nrow(overlap_summary)
