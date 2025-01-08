### Rscript cmplot.R /christen_vcf/vcftools_weir_fst/hybridzone christen G01_G02 G06_G07 G08_G09 ###
### Rscript cmplot.R /verena_vcf/vcftools_weir_fst/hybridzone verena G01_G02 G06_G07 G08_G09 ###
### Rscript cmplot.R /angsd_fst angsd G01_G02 G06_G07 G08_G09 ###

library(openxlsx)
library(CMplot)
library(dplyr)
library(psych)
library(sysfonts) #Rfont
library(showtextdb) #Rfont
library(showtext) #Rfont
par(family  = "Arial") #Rfont
showtext_auto() #Rfont
#conda R unable to find fonts in system. Using this as an alternative, otherwise redundant  

args <- commandArgs(trailingOnly = TRUE)

pat <- c("chr1","chr1A","chr2","chr3","chr4","chr4A","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12",
"chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chr23","chr24","chr25",
"chr26","chr27","chr28","chrZ","chrUN_1","chrUN_2","chrUN_3","chrUN_4","chrUN_5","chrUN_6","chrUN_7","chrUN_8","chrUN_9",
"chrUN_10","chrUN_11","chrUN_12","chrUN_13","chrUN_14","chrUN_15","chrUN_16","chrUN_17")

for (filename in args[3:5]) 
{
cat(filename)
setwd("/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/01_probes/fst")
setwd(paste0(getwd(), args[1]))
x <- read.delim(file=paste(filename,"fst.chr", sep="."), header=TRUE, sep = "")
sorted.x <- x %>% arrange(factor(x$CHROM, levels = pat))
sorted.x.chr <- sorted.x[!grepl("scaffold", sorted.x$CHROM),] #remove scaffolds not mapped to chr
sorted.x.fst <- dfOrder(x,c(-4,2)) #sort dcreasingly by Fst then increasingly by ChrNo.
sorted.x.fst.autosomes <- sorted.x.fst[!grepl("chrZ", sorted.x$CHROM),] # only autosomes
sorted.x.fst.chrz <- sorted.x.fst[grepl("chrZ", sorted.x$CHROM),] # only chrZ
setwd("/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/01_probes/fst/fst_summary")
write.table(sorted.x.fst, file=paste(args[2],filename,"chr_sorted.txt", sep="_"), sep = "\t", row.names = FALSE, col.names = TRUE, quote=FALSE)
write.table(sorted.x.fst.autosomes, file=paste(args[2],filename,"chr_autosomes_only.txt", sep="_"), sep = "\t", row.names = FALSE, col.names = TRUE, quote=FALSE)
write.table(sorted.x.fst.chrz, file=paste(args[2],filename,"chr_chrz_only.txt", sep="_"), sep = "\t", row.names = FALSE, col.names = TRUE, quote=FALSE)
CMplot(sorted.x.chr[,c(1:3,4)],plot.type="m",LOG10=FALSE,col=c("grey30","grey60"),ylim=c(0,1),file="jpg",
memo=paste(filename,args[2],"chr_sorted_only",sep="_"),dpi=150,file.output=TRUE,verbose=TRUE,cex=0.5,width=14,height=6,threshold.lty=1,threshold.lwd=3,threshold.col="red",ylab="Fst")
}
