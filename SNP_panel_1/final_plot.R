### sbatch 4.1_snp_panel_check.sh set6up ###
### Rscript final_plot.R $1 G01_G02 G06_G07 G08_G09 05_final_backup ###
# Rversion 4.0.2 on R studio, R version 3.6.1 on conda

library(devtools)
library(ggplot2)
library(SNPRelate)
library(CMplot)
library(openxlsx)
library(dplyr)
library(psych)
library(ggpubr)
library(sysfonts)
library(showtextdb)
library(showtext)
library(chromPlot) #packageVersion 1.14 on R3.6.1 does not work well

args <- commandArgs(trailingOnly = TRUE)

setwd("/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/01_probes/snp_panel_summary")
setwd(args[5])

# PCA with SNPRelate #
poplist <- read.delim("/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/01_probes/fst/poplist_118.txt") 
colorBlindBlack8  <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")

vcf.fn <- (file=paste("snp_panel",args[1],"final.recode.vcf",sep="_"))
snpgdsVCF2GDS(vcf.fn, "vcf.gds", method="copy.num.of.ref")
snpgdsSummary("vcf.gds")
vcf.gdsfile <- snpgdsOpen("vcf.gds")
vcf.pca <- snpgdsPCA(vcf.gdsfile, autosome.only=FALSE) ##use autosome.only=false to include all loci
names(vcf.pca)
#variance proportion (%)
pc.percent <- vcf.pca$varprop*100
head(round(pc.percent, 2))
print(pc.percent)
tab <- data.frame(sample.id = vcf.pca$sample.id,
EV1 = vcf.pca$eigenvect[,1], # the first eigenvector
EV2 = vcf.pca$eigenvect[,2], # the second eigenvector
EV3 = vcf.pca$eigenvect[,3], 
EV4 = vcf.pca$eigenvect[,4], 
EV5 = vcf.pca$eigenvect[,5], 
EV6 = vcf.pca$eigenvect[,6], 
EV7 = vcf.pca$eigenvect[,7], 
EV8 = vcf.pca$eigenvect[,8], 
EV9 = vcf.pca$eigenvect[,9], 
EV10 = vcf.pca$eigenvect[,10], 
stringsAsFactors = FALSE)
#extract vectors out for external plotting
tab$POP <- poplist$pop
write.csv(tab, file=paste("snp_panel",args[1],"snprelate.txt",sep="_"))

pdf("snp_panel_all_final_pca.pdf")
ggplot(tab, group=POP) +
  geom_point(aes(x=EV1, y=EV2, color=POP, shape=POP), size=4) +
  scale_shape_manual(name="Population", values=c(15, 16, 17, 18, 3, 4, 9, 10, 11, 12, 13, 5)) +
  scale_color_manual(name="Population", values=safe_colorblind_palette) +
  labs(x="PC1", y="PC2") +
  guides(color=guide_legend("Population"),fill=guide_legend("Population")) +
  ggtitle("SNP_panel_all_final") +
  theme_bw(base_size=12) +
  theme(axis.line = element_line(colour = "black"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank(),
  axis.title.x = element_text(size=15),
  axis.title.y = element_text(size=15),
  axis.text.x=element_text (size=12),
  axis.text.y=element_text (size=12),
  legend.text=element_text(size=12),
  legend.title=element_text(size=12)) 
dev.off()

vcf.fn2 <- (file=paste("snp_panel_neutral",args[1],"final.recode.vcf",sep="_"))
snpgdsVCF2GDS(vcf.fn2, "vcf.gds2", method="copy.num.of.ref")
snpgdsSummary("vcf.gds2")
vcf.gdsfile2 <- snpgdsOpen("vcf.gds2")
vcf.pca2 <- snpgdsPCA(vcf.gdsfile2, autosome.only=FALSE) 
names(vcf.pca2)
pc.percent2 <- vcf.pca2$varprop*100
head(round(pc.percent2, 2))
print(pc.percent2)
tab2 <- data.frame(sample.id = vcf.pca2$sample.id,
EV1 = vcf.pca2$eigenvect[,1], # the first eigenvector
EV2 = vcf.pca2$eigenvect[,2], # the second eigenvector
EV3 = vcf.pca2$eigenvect[,3],
EV4 = vcf.pca2$eigenvect[,4],
EV5 = vcf.pca2$eigenvect[,5],
EV6 = vcf.pca2$eigenvect[,6],
EV7 = vcf.pca2$eigenvect[,7],
EV8 = vcf.pca2$eigenvect[,8],
EV9 = vcf.pca2$eigenvect[,9],
EV10 = vcf.pca2$eigenvect[,10],
stringsAsFactors = FALSE)
tab2$POP <- poplist$pop
write.csv(tab2, file=paste("snp_panel_neutral",args[1],"snprelate.txt",sep="_"))

pdf("snp_panel_neutral_final_pca.pdf")
ggplot(tab2, group=POP) +
  geom_point(aes(x=EV1, y=EV2, color=POP, shape=POP), size=4) +
  scale_shape_manual(name="Population", values=c(15, 16, 17, 18, 3, 4, 9, 10, 11, 12, 13, 5)) +  
  scale_color_manual(name="Population", values=safe_colorblind_palette) +
  labs(x="PC1", y="PC2") +
  guides(color=guide_legend("Population"),fill=guide_legend("Population")) +
  ggtitle("SNP_panel_neutral_final") +
  theme_bw(base_size=12) +
  theme(axis.line = element_line(colour = "black"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank(),
  axis.title.x = element_text(size=15),
  axis.title.y = element_text(size=15),
  axis.text.x=element_text (size=12),
  axis.text.y=element_text (size=12),
  legend.text=element_text(size=12),
  legend.title=element_text(size=12)) 
dev.off()

# cmplot #
par(family  = "Arial")
showtext_auto()
pat <- c("chr1","chr1A","chr2","chr3","chr4","chr4A","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chr23","chr24","chr25","chr26","chr27","chr28","chrZ","chrUN_1","chrUN_2","chrUN_3","chrUN_4","chrUN_5","chrUN_6","chrUN_7","chrUN_8","chrUN_9","chrUN_10","chrUN_11","chrUN_12","chrUN_13","chrUN_14","chrUN_15","chrUN_16","chrUN_17")
for (filename in args[2:4]) 
{
cat(filename)
x <- read.delim(file=paste(filename,"weir.fst.chr", sep="."), header=TRUE, sep = "")
sorted.x <- x %>% arrange(factor(x$CHROM, levels = pat))
sorted.x.chr <- sorted.x[!grepl("scaffold", sorted.x$CHROM),] #remove scaffolds not mapped to chr
CMplot(sorted.x.chr[,c(1:3,4)],plot.type="m",LOG10=FALSE,col=c("grey30","grey60"),ylim=c(0,1),file="jpg",
memo=paste(filename,"snp_panel_final_chr_sorted_only",sep="_"),dpi=150,file.output=TRUE,verbose=TRUE,cex=0.5,width=14,height=6,threshold.lty=1,threshold.lwd=3,threshold.col="red",ylab="Fst")
}

# HET plot #
hetbase <- read.delim("../03_neutral/all_neutral_sites_nochrz.het")
het <- read.delim(file=paste("snp_panel_neutral",args[1],"final.het",sep="_"))
corr <- data.frame(all=hetbase$F, snp_panel=het$F)
corr$POP <- poplist$pop
corr <- corr[-c(24), ] #remove anomaly sample A26
pdf("snp_panel_neutral_final_het.pdf")
ggplot(corr) +
geom_point(aes(x=all, y=snp_panel, color=POP, shape=POP), size=2) +
  labs(x = "All_neutral_sites", y = "SNP_panel_neutral_sites") +
  geom_abline(intercept = 0, slope = 1, color="black", linetype="dashed", size=1) +
  scale_color_manual(name="Population", values=safe_colorblind_palette) +
  scale_shape_manual(name="Population", values=c(15, 16, 17, 18, 3, 4, 9, 10, 11, 12, 13, 5)) +
  ggtitle("F_inbreeding_coefficient [1-(FreqObsHet/FreqExHet)]")+
  theme(plot.title = element_text(size=10))
dev.off()

# chrom plot #
# version 1.14 does not allow skipping of missing bands. this will not work by running on conda.
dat1 <- read.delim("snp_panel_outlier_genic_neutral_chrplot.txt", sep="")
dat1[1] <- gsub("scaffold_", "", dat1$Chrom)
###increase size of snps by 10^2 for display purpose
dat1$Start <- dat1$Start-1e2
dat1$End <- dat1$End+1e2
###colours
dat1$Colors <-
  ifelse(dat1$Group == "dxy","black",
  ifelse(dat1$Group == "dxy_genic","grey1",
  ifelse(dat1$Group == "dxy_knief","grey2",
  ifelse(dat1$Group == "fst", "grey3",
  ifelse(dat1$Group == "fst_dxy", "grey4",
  ifelse(dat1$Group == "fst_dxy_genic", "grey5", 
  ifelse(dat1$Group == "fst_dxy_knief", "grey6",
  ifelse(dat1$Group == "fst_dxy_knief_genic", "grey7",
  ifelse(dat1$Group == "fst_genic", "grey8",
  ifelse(dat1$Group == "fst_knief", "grey9",
  ifelse(dat1$Group == "fst_knief_genic", "grey10",
  ifelse(dat1$Group == "genic",  "orange",
  ifelse(dat1$Group == "knief","red1",
  ifelse(dat1$Group == "knief_genic","red2",
  ifelse(dat1$Group == "neutral", "royalblue","royalblue")))))))))))))))
pdf('chrplot_1e2_1to50.pdf')
chromPlot(bands=dat1, chr=c(0:50))
dev.off()
pdf('chrplot_1e2_51to100.pdf')
chromPlot(bands=dat1, chr=c("M",51:100))
dev.off()
#grouping by chromosomes mapped to ref5.5
chr18 <- dat1[(dat1$ChrNo=="chr18"),]
chr21 <- dat1[(dat1$ChrNo=="chr21"),]
chr23 <- dat1[(dat1$ChrNo=="chr23"),]
pdf('chrplot_1e2_chr18.pdf')
chromPlot(bands=chr18)
dev.off()
pdf('chrplot_1e2_chr21.pdf')
chromPlot(bands=chr21)
dev.off()
pdf('chrplot_1e2_chr23.pdf')
chromPlot(bands=chr23)
dev.off()

