library(ggplot2)
library(readxl)
library(dplyr)
library(gridExtra)
library(grid)
library(tibble)
library(VennDiagram)
library(data.table)
library(tidyr)
library(ggsignif)

setwd("/Users/chyiyin/Dropbox/CORVID_baits/Analyses")
dat <- read.delim("baitscomparison23.txt",header=TRUE, sep="\t") 

#data grouped by country and sorted from the youngest to oldest samples
dat$perdedupmappedreads <- (dat$Nr..Dedup..Mapped.Reads - dat$Reads.to.be.ommitted) / (dat$Nr..Reads.Into.Mapping - dat$Reads.to.be.ommitted) *100 
dat$complexity <- dat$Mean.cov.of.target.104k.SNPsite / dat$Expected.genomic.coverage..based.on.qPCR.

dat1 = dat %>% group_by(Country) %>% arrange(desc(Age), .by_group = TRUE)
twistdat = dat1 %>% filter(Method == "Twist") 
mybaitsdat = dat1 %>% filter(Method == "myBaits")

##### Figure 2 scatterplot #####
par(mfcol=c(2,4))

#1
x=twistdat$perdedupmappedreads
y=mybaitsdat$perdedupmappedreads
plot(x=twistdat$perdedupmappedreads, y=mybaitsdat$perdedupmappedreads,main = "Deduplicate reads mapped to\n crow genome (%)",
     xlab = "Twist", ylab = "myBaits", pch = 19, frame = FALSE,xlim=c(0,50), ylim=c(0,50), cex.lab = 1.4, cex.axis=1.3)
abline(lm(y ~ x), col = "blue")
abline(0,1,col = "black",lty=2)
mtext("A", side = 3, adj = 0, line = 1.5, cex = 1.2, font = 2)

#2
x=twistdat$complexity
y=mybaitsdat$complexity
plot(x=twistdat$complexity, y=mybaitsdat$complexity,main = "Complexity\n (obs vs exp coverage)",
     xlab = "Twist", ylab = "myBaits", pch = 19, frame = FALSE, xlim=c(0,1), ylim=c(0,1), cex.lab = 1.4, cex.axis=1.3)
abline(lm(y ~ x), col = "blue")
abline(0,1,col = "black",lty=2)
mtext("E", side = 3, adj = 0, line = 1.5, cex = 1.2, font = 2)

#3
x=twistdat$On.target.alignment.adjusted.to.104K.80bp.panel....
y=mybaitsdat$On.target.alignment.adjusted.to.104K.80bp.panel....
plot(x=twistdat$On.target.alignment.adjusted.to.104K.80bp.panel...., y=mybaitsdat$On.target.alignment.adjusted.to.104K.80bp.panel....,main = "On-target deduplicate alignment to\n 104K_80bp panel (%)",
     xlab = "Twist", ylab = "myBaits", pch = 19, frame = FALSE, xlim=c(0,70), ylim=c(0,70), cex.lab = 1.4, cex.axis=1.3)
abline(lm(y ~ x), col = "blue")
abline(0,1,col = "black",lty=2)
mtext("B", side = 3, adj = 0, line = 1.5, cex = 1.2, font = 2)

#4
x=twistdat$MT.to.Nuclear.Ratio
y=mybaitsdat$MT.to.Nuclear.Ratio
plot(x=twistdat$MT.to.Nuclear.Ratio, y=mybaitsdat$MT.to.Nuclear.Ratio,main = "MT to Nuclear Ratio\n (coverage)",
     xlab = "Twist", ylab = "myBaits", pch = 19, frame = FALSE, xlim=c(0,100), ylim=c(0,100), cex.lab = 1.4, cex.axis=1.3)
abline(lm(y ~ x), col = "blue")
abline(0,1,col = "black",lty=2)
mtext("F", side = 3, adj = 0, line = 1.5, cex = 1.2, font = 2)

#5
x=twistdat$On.target.alignment.to.original.bait.panel....
y=mybaitsdat$On.target.alignment.to.original.bait.panel....
plot(x=twistdat$On.target.alignment.to.original.bait.panel...., y=mybaitsdat$On.target.alignment.to.original.bait.panel....,main = "On-target deduplicate alignment to\n original probe length & size (%)",
     xlab = "Twist", ylab = "myBaits", pch = 19, frame = FALSE, xlim=c(0,70), ylim=c(0,70), cex.lab = 1.4, cex.axis=1.3)
abline(lm(y ~ x), col = "blue")
abline(0,1,col = "black",lty=2)
mtext("C", side = 3, adj = 0, line = 1.5, cex = 1.2, font = 2)

#6
x=twistdat$Mean.Length.Mapped.Reads
y=mybaitsdat$Mean.Length.Mapped.Reads
plot(x=twistdat$Mean.Length.Mapped.Reads, y=mybaitsdat$Mean.Length.Mapped.Reads,main = "Mean length of deduplicate\n mapped reads (bp)",
     xlab = "Twist", ylab = "myBaits", pch = 19, frame = FALSE, xlim=c(0,70), ylim=c(0,70), cex.lab = 1.4, cex.axis=1.3)
abline(lm(y ~ x), col = "blue")
abline(0,1,col = "black",lty=2)
mtext("G", side = 3, adj = 0, line = 1.5, cex = 1.2, font = 2)

#7
x=twistdat$Mean.cov.of.target.104k.SNPsite
y=mybaitsdat$Mean.cov.of.target.104k.SNPsite
plot(x=twistdat$Mean.cov.of.target.104k.SNPsite, y=mybaitsdat$Mean.cov.of.target.104k.SNPsite,main = "Mean coverage of\n 104K SNP sites",
     xlab = "Twist", ylab = "myBaits", pch = 19, frame = FALSE, xlim=c(0,30), ylim=c(0,30), cex.lab = 1.4, cex.axis=1.3)
abline(lm(y ~ x), col = "blue")
abline(0,1,col = "black",lty=2)
mtext("D", side = 3, adj = 0, line = 1.5, cex = 1.2, font = 2)

#8
x=twistdat$X5.Prime.C.T.1st.base.on.target.104k_80bpdedup
y=mybaitsdat$X5.Prime.C.T.1st.base.on.target.104k_80bpdedup
plot(x=twistdat$X5.Prime.C.T.1st.base.on.target.104k_80bpdedup, y=mybaitsdat$X5.Prime.C.T.1st.base.on.target.104k_80bpdedup,main = "C-T subst on 1st bp of\n on-target reads (%)",
     xlab = "Twist", ylab = "myBaits", pch = 19, frame = FALSE, xlim=c(0,60), ylim=c(0,60), cex.lab = 1.4, cex.axis=1.3)
abline(lm(y ~ x), col = "blue")
abline(0,1,col = "black",lty=2)
mtext("H", side = 3, adj = 0, line = 1.5, cex = 1.2, font = 2)

dev.off() #saved 14x7

##### Figure 3 Boxplot #####
count_80foldcov1 <- dat1 %>%
  filter( at.least.1X.of.target.104k.SNPsite.... > 79) %>% 
  group_by(Method) %>% 
  summarize(count_above_80 = n())
count_80foldcov2 <- dat1 %>%
  filter( at.least.2X.of.target.104k.SNPsite.... > 79) %>% 
  group_by(Method) %>% 
  summarize(count_above_80 = n())
count_80foldcov3 <- dat1 %>%
  filter( at.least.3X.of.target.104k.SNPsite.... > 79) %>% 
  group_by(Method) %>% 
  summarize(count_above_80 = n())
count_80foldcov4 <- dat1 %>%
  filter( at.least.4X.of.target.104k.SNPsite.... > 79) %>% 
  group_by(Method) %>% 
  summarize(count_above_80 = n())

p1<- ggplot(dat1, aes(x = Method, y =  at.least.1X.of.target.104k.SNPsite...., fill = Method)) +
  geom_violin(trim = FALSE, alpha = 0.5, show.legend = FALSE) +
  geom_boxplot(width=0.05, alpha = 0.1, show.legend = FALSE) +
  geom_signif(comparisons = list(c("myBaits", "Twist")),map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05),tip_length = 0,y_position=0.7) +
  geom_point(position = position_jitter(width = 0.2), size = 0.8, alpha = 0.6, show.legend = FALSE) +
  labs(title="At least 1X coverage",y = "Percentage of target SNPs (%)", x="") +
  coord_cartesian(ylim = c(0, 100)) +
  geom_text(data = count_80foldcov1, aes(x=Method,y = 80, label = c(paste("N = ",count_80foldcov1[1,2],sep=""), paste("N = ",count_80foldcov1[2,2],sep=""))), 
            vjust = -0.5, hjust=1.2, size = 4, color = "black") +
  geom_hline(yintercept = 80, linetype = "dashed", color = "black") +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 10)) + 
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

p2<- ggplot(dat1, aes(x = Method, y =  at.least.2X.of.target.104k.SNPsite...., fill = Method)) +
  geom_violin(trim = FALSE, alpha = 0.5, show.legend = FALSE) +
  geom_boxplot(width=0.05, alpha = 0.1, show.legend = FALSE) +
  geom_signif(comparisons = list(c("myBaits", "Twist")),map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05),tip_length = 0,y_position=0.7) +
  geom_point(position = position_jitter(width = 0.2), size = 0.8, alpha = 0.6, show.legend = FALSE) + 
  labs(title="At least 2X coverage",y = "", x="") +
  geom_text(data = count_80foldcov2, aes(x=Method,y = 80, label = c(paste("N = ",count_80foldcov2[1,2],sep=""), paste("N = ",count_80foldcov2[2,2],sep=""))), 
            vjust = -0.5, hjust=1.2, size = 4, color = "black") +
  geom_hline(yintercept = 80, linetype = "dashed", color = "black") +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 10)) + 
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

p3<- ggplot(dat1, aes(x = Method, y =  at.least.3X.of.target.104k.SNPsite...., fill = Method)) +
  geom_violin(trim = FALSE, alpha = 0.5, show.legend = FALSE) +
  geom_boxplot(width=0.05, alpha = 0.1, show.legend = FALSE) +
  geom_signif(comparisons = list(c("myBaits", "Twist")),map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05),tip_length = 0,y_position=0.7) +
  geom_point(position = position_jitter(width = 0.2), size = 0.8, alpha = 0.6, show.legend = FALSE) +  # Add jittered points
  labs(title="At least 3X coverage",y = "Percentage of target SNPs (%)", x="") +
  geom_text(data = count_80foldcov3, aes(x=Method,y = 80, label = c(paste("N = ",count_80foldcov3[1,2],sep=""), paste("N = ",count_80foldcov3[2,2],sep=""))), 
            vjust = -0.5, hjust=1.2, size = 4, color = "black") +
  geom_hline(yintercept = 80, linetype = "dashed", color = "black") +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 10)) + 
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

p4<- ggplot(dat1, aes(x = Method, y =  at.least.4X.of.target.104k.SNPsite...., fill = Method)) +
  geom_violin(trim = FALSE, alpha = 0.5, show.legend = FALSE) +
  geom_boxplot(width=0.05, alpha = 0.1, show.legend = FALSE) +
  geom_signif(comparisons = list(c("myBaits", "Twist")),map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05),tip_length = 0,y_position=0.7) +
  geom_point(position = position_jitter(width = 0.2), size = 0.8, alpha = 0.6, show.legend = FALSE) +  # Add jittered points
  labs(title="At least 4X coverage",y = "", x="") +
  geom_text(data = count_80foldcov4, aes(x=Method,y = 80, label = c(paste("N = ",count_80foldcov4[1,2],sep=""), paste("N = ",count_80foldcov4[2,2],sep=""))), 
            vjust = -0.5, hjust=1.2, size = 4, color = "black") +
  geom_hline(yintercept = 80, linetype = "dashed", color = "black") +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 10)) + 
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

grid.arrange(p1,p2,p3,p4, nrow=2) #saved 6x6

##### Supplementary figures SNP coverage distribution from samtools depth #####
library(purrr)
mybaitscov <- read.delim("coverage_104k_mybaits.txt",header=FALSE, sep="\t") 
twistcov <- read.delim("coverage_104k_twist.txt",header=FALSE, sep="\t") 
popbaits <- read.table("/Users/chyiyin/Dropbox/CORVID_baits/Analyses/angsd/popbaits.txt", sep="\t", header=TRUE) 
popbaits_selected <- popbaits %>% 
  filter(Country %in% c("PL","B")) %>%
  filter(Samples != "DVT016" & Samples != "DVT022" & Samples != "KZR002")

mybaitscov <- mybaitscov %>%
  set_names(c("scaffold", "pos", popbaits$Samples))
twistcov <- twistcov %>%
  set_names(c("scaffold", "pos", popbaits$Samples))

plot_list <- list()
for (type in c("mybaitscov", "twistcov")) {
  for (i in popbaits$Samples) {
    if (i %in% names(get(type))) {
      fill_color <- ifelse(type == "mybaitscov", "#F8766D", "#00BFC4")
      plot <- ggplot(data=get(type), aes_string(x = i)) +
        geom_histogram(binwidth=1, fill = fill_color, color = "black") +
        labs(title = i, x = "Coverage", y = "No. of SNP sites") +
        scale_x_continuous(limits = c(0, 20)) +
        theme_minimal() + 
        theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
              panel.grid.minor = element_blank(),
              axis.title.x = element_text(size = 8, face = "bold"),
              axis.title.y = element_text(size = 8, face = "bold"),
              axis.text.x = element_text(size = 8, face = "bold"),
              axis.text.y = element_text(size = 8, face = "bold"))
      
      plot_list[[i]] <- plot
    } else {
      warning(paste("Column", i, "not found in mybaitscov"))
    } } 
  pdf(file=paste("coverage_distribution_104k_",type,"_selected",sep=""))
  grid.arrange(grobs=plot_list[popbaits_selected$Samples], ncol=4)
  dev.off()
}

#### Figure 4 normalized coverage across GC content ####
mybaitsgc <- read.delim("./gc/GC_coverage_summary_mybaits.txt",header=TRUE, sep="\t") 
twistgc <- read.delim("./gc/GC_coverage_summary_twist.txt",header=TRUE, sep="\t") 
names(mybaitsgc) <- gsub("_mybaits", "", names(mybaitsgc))
names(twistgc) <- gsub("_TE", "", names(twistgc))
mybaitsgc_selected <- mybaitsgc %>% 
  select(c("GC", popbaits_selected$Samples)) %>%
  pivot_longer(cols = popbaits_selected$Samples, names_to = "variable", values_to = "myBaits") %>%
  filter(GC>0)
twistgc_selected <- twistgc %>% 
  select(c("GC", popbaits_selected$Samples)) %>%
  pivot_longer(cols = popbaits_selected$Samples, names_to = "variable", values_to = "Twist") %>%
  filter(GC>0)
merged_gc <- cbind(mybaitsgc_selected,Twist=twistgc_selected$Twist) %>%
  pivot_longer(cols = c(myBaits, Twist), names_to = "Method", values_to = "y")

gccov <- ggplot(merged_gc, aes(x = GC, y = y, color = Method)) +
  geom_point(size = 0.5,alpha=0.4) +  
  geom_smooth(method = "gam", se = TRUE) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 5)) +
  scale_color_manual(values = c("#F8766D", "#00BFC4"), name = "Method", labels = c("myBaits", "Twist")) +
  labs(title = "",x = "GC content",y = "Normalized coverage") +
  theme_minimal() +
  theme(legend.position = c(0.9, 0.8),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))

#### Figure 4 Boxplot for GC / AT dropout ####
dropout <- read.delim("./gc/GC_AT_dropout_summary.txt",header=TRUE, sep="\t") 
dropout_selected <- dropout %>% separate(col=Sample,into=c("Sample","Method"),sep="_") %>%  
  mutate(Method=recode(Method, "TE"="Twist", "mybaits"="myBaits")) %>%
  filter(Sample %in% popbaits_selected$Samples) 

drop1<- ggplot(dropout_selected, aes(x = Method, y = GC_Dropout, fill = Method)) +
  geom_violin(trim = FALSE, alpha = 0.5, show.legend = FALSE) +
  geom_boxplot(width=0.05, alpha = 0.1, show.legend = FALSE) +
  geom_signif(comparisons = list(c("myBaits", "Twist")),map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05),tip_length = 0,y_position=1) +
  geom_point(position = position_jitter(width = 0.2), size = 0.8, alpha = 0.6, show.legend = FALSE) + 
  labs(title="GC dropout",y = "", x="") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
drop2<- ggplot(dropout_selected, aes(x = Method, y = AT_Dropout, fill = Method)) +
  geom_violin(trim = FALSE, alpha = 0.5, show.legend = FALSE) +
  geom_boxplot(width=0.05, alpha = 0.1, show.legend = FALSE) +
  geom_signif(comparisons = list(c("myBaits", "Twist")),tip_length = 0, y_position=25,map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05)) +
  geom_point(position = position_jitter(width = 0.2), size = 0.8, alpha = 0.6, show.legend = FALSE) + 
  labs(title="AT dropout",y = "", x="") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
grid.arrange(gccov, arrangeGrob(drop2,drop1, ncol=2)) #saved 8x8

##### Figure 5 SNP quality #####
dat1 = dat %>% group_by(Country) %>% arrange(desc(Age), .by_group = FALSE) %>% 
  filter(Country %in% c("PL","B")) %>%
  filter(Sample.Name != "DVT016" & Sample.Name != "DVT022" & Sample.Name != "KZR002")
twistdat = dat1 %>% filter(Method == "Twist") 
mybaitsdat = dat1 %>% filter(Method == "myBaits")
scale_factor <- max(dat1$X5.Prime.C.T.1st.base.on.target.104k_80bpdedup) / max(dat1$Age,na.rm=TRUE)
p5 <- ggplot(data=dat1,aes(x=Sample.Name,y=X5.Prime.C.T.1st.base.on.target.104k_80bpdedup, fill=Method)) + 
  geom_bar(stat="identity",position="dodge") + 
  geom_line(aes(x=Sample.Name,y=Age* scale_factor, group = 1), color = "black", size = 0.8, linetype="dashed") +
  scale_y_continuous(name = "C-T subst on 1st bp of on-target reads (%)", sec.axis = sec_axis(~ . / scale_factor, name = "Age (years ago)")) +  # Secondary axis for the rate) +
  theme_bw() +
  theme(axis.text.x = element_text(face="bold", angle=90)) +
  scale_x_discrete(limits=mybaitsdat$Sample.Name) +
  theme(axis.title.x=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.9, 0.7),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))

##### SNP mismatch #####
setwd("/Users/chyiyin/Dropbox/CORVID_baits/Analyses/angsd") 
twistsnp <- read.table(gzfile("twist_104kpanel_hapconsensus_maxmis_q20.haplo.gz"), sep = "\t", header=TRUE) #76303
mybaitssnp <- read.table(gzfile("mybaits_104kpanel_hapconsensus_maxmis_q20.haplo.gz"), sep = "\t", header=TRUE) #73038
shotgunsnp <- read.table(gzfile("shotgun_noudg_104kpanel_hapconsensus_maxmis_q20.haplo.gz"), sep = "\t", header=TRUE) #42353
shotgungeno <- read.table(gzfile("shotgun_noudg_104kpanel_geno_maxmis_q20.geno.gz"), sep = "\t", header=FALSE)
shotgungeno <- shotgungeno[-15]
# the filter for snpcalling on angsd is -minMapQ 20 -minQ 20

popbaits <- read.table("popbaits.txt", sep="\t", header=TRUE) 
popbaits_selected <- popbaits %>% 
  filter(Country %in% c("PL","B")) %>%
  filter(Samples != "DVT016" & Samples != "DVT022" & Samples != "KZR002")
popshotgun <- read.table("popshotgun.txt", sep="\t", header=TRUE) 
popshotgun_selected <- popshotgun %>%
  filter(Country %in% c("PL","B")) %>%
  filter(Samples != "DVT016" & Samples != "DVT022" & Samples != "KZR002")

colnames(twistsnp) <-  c("chr","pos","major",popbaits$Samples)
colnames(mybaitssnp) <-  c("chr","pos","major",popbaits$Samples)
colnames(shotgunsnp) <-  c("chr","pos","major",popshotgun$Samples)
colnames(shotgungeno) <-  c("chr","pos",popshotgun$Samples)

twistsnp <- twistsnp %>% mutate(scaff_name=paste(chr,pos,sep="_")) %>% 
  column_to_rownames(var = "scaff_name") %>%
  mutate(across(everything(), ~ ifelse(. == "N", NA, .))) %>%
  select(popbaits_selected$Samples)
mybaitssnp <- mybaitssnp %>% mutate(scaff_name=paste(chr,pos,sep="_"))  %>% 
  column_to_rownames(var = "scaff_name") %>%
  mutate(across(everything(), ~ ifelse(. == "N", NA, .))) %>%
  select(popbaits_selected$Samples)
shotgunsnp <- shotgunsnp %>% mutate(scaff_name=paste(chr,pos,sep="_"))  %>% 
  column_to_rownames(var = "scaff_name") %>%
  mutate(across(everything(), ~ ifelse(. == "N", NA, .))) %>%
  select(popshotgun_selected$Samples)
shotgungeno <- shotgungeno %>% mutate(scaff_name=paste(chr,pos,sep="_"))  %>% 
  column_to_rownames(var = "scaff_name") %>%
  mutate(across(everything(), ~ ifelse(. == "NN", NA, .))) %>%
  select(popshotgun_selected$Samples)

# Venn diagram
set1 <- rownames(mybaitssnp)
set2 <- rownames(twistsnp)
set3 <- rownames(shotgunsnp)
venn.plot <- venn.diagram(x = list("myBaits"=set1, "Twist"=set2, "shotgun"=set3),
                          category.names = c("myBaits", "Twist", "shotgun"),
                          filename = NULL,  # To display in R instead of saving as a file
                          margin = 0.05,
                          output = TRUE,
                          fill = c("#F8766D", "#00BFC4", "orange"),
                          alpha = 0.5,
                          cex = 1,
                          cat.fontfamily = "sans",
                          fontfamily = "sans",
                          fontface = "bold",
                          cat.fontface = "bold")


# Subset both files to only common rows and other filters
common_rows <- intersect(rownames(mybaitssnp), rownames(twistsnp)) #67681
twistsnp_common <- twistsnp[common_rows, ] 
mybaitssnp_common <- mybaitssnp[common_rows, ] 
if (!all(dim(twistsnp_common) == dim(mybaitssnp_common))) {
  stop("The files do not have the same number of common rows and columns.")
}

# Compare values in common rows
differences <- twistsnp_common != mybaitssnp_common  # This creates a logical matrix of TRUE/FALSE
diff_rows <- rowSums(differences, na.rm = TRUE) > 0  # Identify rows with at least one difference
mismatch_counts <- sort(colSums(differences, na.rm = TRUE), decreasing=TRUE)

# Check againts shotgun SNP
mismatchresults <- data.frame(Sample=popshotgun_selected$Samples,myBaits_match_shotgun=1,Twist_match_shotgun=2,heterozygosity=3)
for (i in 1:length(popshotgun_selected$Samples)) {
  sample_name <- popshotgun_selected$Samples[i] 
  mismatched_rows <- which(differences[,sample_name], arr.ind = TRUE)
  mismatched_row_names <- rownames(differences)[mismatched_rows]
  
  # Check against shotgun
  check_shotgun <- intersect(mismatched_row_names, rownames(shotgunsnp))
  twist2 <- twistsnp_common[check_shotgun,] %>%
    select(sample_name)
  mybaits2 <- mybaitssnp_common[check_shotgun,] %>%
    select(sample_name)
  shotgun2 <- shotgunsnp[check_shotgun,] %>%
    select(sample_name)
  geno <- shotgungeno[check_shotgun,] %>%
    select(sample_name)
  twist2 <- tibble::rownames_to_column(twist2, var = "rowname")
  mybaits2 <- tibble::rownames_to_column(mybaits2, var = "rowname")
  shotgun2 <- tibble::rownames_to_column(shotgun2, var = "rowname")
  geno <- tibble::rownames_to_column(geno, var = "rowname")
  geno2 <- geno %>% separate(sample_name, into =c("allele1","allele2"),sep=1) %>%
    mutate(het=(allele1!=allele2)*1)
  het <- sum(geno2$het, na.rm = TRUE)/sum(!is.na(geno2$het))*100
  
  merged_df <- mybaits2 %>%
    full_join(twist2, by = "rowname") %>%
    full_join(shotgun2, by = "rowname")
  colnames(merged_df) <- c("rowname","myBaits","Twist","shotgun")
  
  mybaits_match_shotgun <- sum(merged_df$myBaits == merged_df$shotgun, na.rm=TRUE) 
  mismatchresults[i,2] <- mybaits_match_shotgun
  twist_match_shotgun <- sum(merged_df$Twist == merged_df$shotgun, na.rm=TRUE) 
  mismatchresults[i,3] <- twist_match_shotgun
  mismatchresults[i,4] <- het
}

mismatchresults2 <- melt(as.data.table(mismatchresults), id.vars=c(1,4), variable.name="Method", value.name="Matches to shotgun")
m1 <- ggplot(data=mismatchresults2,aes(x=Sample,y=`Matches to shotgun`, fill=Method)) + 
  geom_bar(stat="identity",position="dodge",show.legend = FALSE) + 
  geom_text(aes(label = `Matches to shotgun`),position = position_dodge(width = 0.9), vjust = -0.3, size=3) +
  geom_text(aes(y = -max(`Matches to shotgun`)*0.05,label = paste("het=",format(heterozygosity,digits=3,scientific=FALSE),"%",sep="")),size=3.5) +
  theme_bw() +
  scale_fill_discrete(labels = c("myBaits", "Twist")) +
  theme(axis.text.x = element_text(face="bold")) +
  theme(axis.title.x=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 

grid.arrange(p5,arrangeGrob(venn.plot, m1, ncol = 2))

