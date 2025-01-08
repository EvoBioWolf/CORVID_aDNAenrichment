### Rscript outlier_genic_neutral $10 ###
### $10 = set6up ###

library(dplyr)
library(devtools)
library(grid)
library(ggplot2)
#library(ggvenn) #only on Rstudio
library(VennDiagram)
library(base)

setwd("/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/01_probes/snp_panel_summary/03_neutral")
args <- commandArgs(trailingOnly = TRUE)

outlier <- read.delim("../01_outlier/all_chrz_fst_dxy99.95_knief_strict_25k_biallelic_only.txt", header=TRUE, sep="")
genic <- read.delim("../02_genic/group1genes_G01G02_14k_biallelic_only.txt", header=TRUE, sep="")
neutral <- read.delim(file=paste(args[1],"biallelic_reduced_positions.txt",sep="_"), header=TRUE, sep="")
outlier_genic_neutral <- list(outlier = as.character(outlier$LocusName), genic = as.character(genic$LocusName), neutral = as.character(neutral$LocusName))

venn.plot <- venn.diagram(outlier_genic_neutral,
filename = NULL,
fill = c("#0073C2FF", "#EFC000FF", "#868686FF"),
main="outlier_genic_neutral")
pdf(file=paste(args[1],"outlier_genic_neutral.pdf",sep="_"))
grid.draw(venn.plot)
dev.off()

overlap <- calculate.overlap(outlier_genic_neutral)
summary(overlap)
overlap_a5 <- data.frame(LocusName=overlap$a5, OverlapType="outlier_genic_neutral")
overlap_a2 <- data.frame(LocusName=overlap$a2, OverlapType="outlier_genic")
overlap_a4 <- data.frame(LocusName=overlap$a4, OverlapType="outlier") #treating outlier_neutral sites as outliers
overlap_a6 <- data.frame(LocusName=overlap$a6, OverlapType="genic_neutral")
overlap_a1 <- data.frame(LocusName=overlap$a1, OverlapType="outlier")
overlap_a3 <- data.frame(LocusName=overlap$a3, OverlapType="genic")
overlap_a7 <- data.frame(LocusName=overlap$a7, OverlapType="neutral")
overlap_summary <- merge(overlap_a5, merge(overlap_a2, merge(overlap_a4, merge(overlap_a6, merge(overlap_a1,
merge(overlap_a3, overlap_a7, all=TRUE), all=TRUE), all=TRUE), all=TRUE), all=TRUE), all=TRUE)
overlap_summary_all <- as.data.frame(do.call(rbind, strsplit(as.character(overlap_summary$LocusName), split=":")))
overlap_summary$Scaffold <- overlap_summary_all$V1
overlap_summary$Pos <- overlap_summary_all$V2
write.table(overlap_summary, file=paste(args[1],"outlier_genic_neutral.txt",sep="_"), sep="\t", row.names=FALSE, col.names=TRUE, 
quote=FALSE)
nrow(overlap_summary)
