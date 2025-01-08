library(dplyr)
library(devtools)
library(grid)
library(ggplot2)
#library(ggvenn) #only available on Rstudio
library(VennDiagram)
library(base)

setwd("/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/01_probes/snp_panel_summary")

dxy99.95 <- read.delim("./01_outlier/G01_G02_persite_dxy99.95.txt", header=TRUE, sep="")
fst99.9filtered <- read.delim("./01_outlier/allhybridzones_persite_chrzinc_fst99.9_filtered.txt", header=TRUE, sep="")
kniefprobes <- read.delim("../goldengate_snps/data_TableS1.txt", header=TRUE, sep="\t")
genic <- read.delim("./02_genic/group1genes_G01G02_14k.txt", header=TRUE, sep="")

outlier_genic <- list(fst99.9filtered = as.character(fst99.9filtered$LocusName), dxy99.95 = as.character(dxy99.95$LocusName),
Kniefprobes = as.character(kniefprobes$LocusName), genic = as.character(genic$LocusName))

### using ggvenn ###
# pdf("outlier_genic.pdf")
# ggvenn(outlier_genic,
# fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
# stroke_size = 0.5, set_name_size = 3) +
# labs(title="fst99.9filtered_dxy99.95_Knief_genic, n=45593") + theme(plot.title = element_text(hjust=0.5))
# dev.off()
###

venn.plot <- venn.diagram(outlier_genic,
filename = NULL,
fill = c("#0073C2FF", "#EFC000FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
main="fst99.9filtered_dxy99.95_Knief_genic")
pdf("outlier_genic.pdf")
grid.draw(venn.plot)
dev.off()

overlap <- calculate.overlap(outlier_genic)
summary(overlap)
overlap_a6 <- data.frame(LocusName=overlap$a6, OverlapType="fst_dxy_knief_genic")
overlap_a12 <- data.frame(LocusName=overlap$a12, OverlapType="fst_dxy_knief")
overlap_a11 <- data.frame(LocusName=overlap$a11, OverlapType="fst_dxy_genic")
overlap_a5 <- data.frame(LocusName=overlap$a5, OverlapType="fst_knief_genic")
overlap_a15 <- data.frame(LocusName=overlap$a15, OverlapType="fst_dxy")
overlap_a4 <- data.frame(LocusName=overlap$a4, OverlapType="fst_knief")
overlap_a10 <- data.frame(LocusName=overlap$a10, OverlapType="fst_genic")
overlap_a13 <- data.frame(LocusName=overlap$a13, OverlapType="dxy_knief")
overlap_a8 <- data.frame(LocusName=overlap$a8, OverlapType="dxy_genic")
overlap_a2 <- data.frame(LocusName=overlap$a2, OverlapType="knief_genic")
overlap_a9 <- data.frame(LocusName=overlap$a9, OverlapType="fst")
overlap_a14 <- data.frame(LocusName=overlap$a14, OverlapType="dxy")
overlap_a1 <- data.frame(LocusName=overlap$a1, OverlapType="knief")
overlap_a3 <- data.frame(LocusName=overlap$a3, OverlapType="genic")
overlap_summary <- merge(overlap_a6, merge(overlap_a12, merge(overlap_a11, merge(overlap_a5, merge(overlap_a15,
merge(overlap_a4, merge(overlap_a10, merge(overlap_a13, merge(overlap_a8, merge(overlap_a2,
merge(overlap_a9, merge(overlap_a14, merge(overlap_a1, overlap_a3, all=TRUE),
all=TRUE), all=TRUE), all=TRUE), all=TRUE), all=TRUE), all=TRUE), all=TRUE), all=TRUE), all=TRUE), all=TRUE), all=TRUE), all=TRUE)
overlap_summary_all <- as.data.frame(do.call(rbind, strsplit(overlap_summary$LocusName, split=":")))
overlap_summary$Scaffold <- overlap_summary_all$V1
overlap_summary$Pos <- overlap_summary_all$V2
write.table(overlap_summary, file="./02_genic/all_chrz_fst_dxy99.95_knief_genic_46k.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
nrow(overlap_summary)
