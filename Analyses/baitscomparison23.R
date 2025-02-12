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

#2 E
x=twistdat$Mean.cov.of.target.104k.SNPsite
y=mybaitsdat$Mean.cov.of.target.104k.SNPsite
plot(x=twistdat$Mean.cov.of.target.104k.SNPsite, y=mybaitsdat$Mean.cov.of.target.104k.SNPsite,main = "Mean coverage of\n 104K SNP sites",
     xlab = "Twist", ylab = "myBaits", pch = 19, frame = FALSE, xlim=c(0,50), ylim=c(0,50), cex.lab = 1.4, cex.axis=1.3)
abline(lm(y ~ x), col = "blue")
abline(0,1,col = "black",lty=2)
mtext("E", side = 3, adj = 0, line = 1.5, cex = 1.2, font = 2)
points_to_label <- c(7,8,11)
text(twistdat$Mean.cov.of.target.104k.SNPsite[points_to_label],mybaitsdat$Mean.cov.of.target.104k.SNPsite[points_to_label], 
     labels = twistdat$Sample.Name[points_to_label], 
     pos = 4, cex = 1, col = "red")

#3
x=twistdat$On.target.alignment.adjusted.to.104K.80bp.panel....
y=mybaitsdat$On.target.alignment.adjusted.to.104K.80bp.panel....
plot(x=twistdat$On.target.alignment.adjusted.to.104K.80bp.panel...., y=mybaitsdat$On.target.alignment.adjusted.to.104K.80bp.panel....,main = "On-target deduplicate alignment to\n 104K_80bp panel (%)",
     xlab = "Twist", ylab = "myBaits", pch = 19, frame = FALSE, xlim=c(0,70), ylim=c(0,70), cex.lab = 1.4, cex.axis=1.3)
abline(lm(y ~ x), col = "blue")
abline(0,1,col = "black",lty=2)
mtext("B", side = 3, adj = 0, line = 1.5, cex = 1.2, font = 2)

#4 F
x=twistdat$complexity
y=mybaitsdat$complexity
plot(x=twistdat$complexity, y=mybaitsdat$complexity,main = "Observed vs expected\n coverage",
     xlab = "Twist", ylab = "myBaits", pch = 19, frame = FALSE, xlim=c(0,1), ylim=c(0,1), cex.lab = 1.4, cex.axis=1.3)
abline(lm(y ~ x), col = "blue")
abline(0,1,col = "black",lty=2)
mtext("F", side = 3, adj = 0, line = 1.5, cex = 1.2, font = 2)
text(twistdat$complexity[points_to_label],mybaitsdat$complexity[points_to_label], 
     labels = twistdat$Sample.Name[points_to_label], 
     pos = 4, cex = 1, col = "red")


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

#7 D
x=twistdat$MT.to.Nuclear.Ratio
y=mybaitsdat$MT.to.Nuclear.Ratio
plot(x=twistdat$MT.to.Nuclear.Ratio, y=mybaitsdat$MT.to.Nuclear.Ratio,main = "MT to Nuclear Ratio\n (coverage)",
     xlab = "Twist", ylab = "myBaits", pch = 19, frame = FALSE, xlim=c(0,100), ylim=c(0,100), cex.lab = 1.4, cex.axis=1.3)
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
shotgungeno <- read.table(gzfile("shotgun_noudg_104kpanel_geno_maxmis_q20_dp3.geno.gz"), sep = "\t", header=FALSE) #104K with Ns
shotgungeno <- shotgungeno[-15]
# the filter for snpcalling on angsd is -minMapQ 20 -minQ 20 and also dp>=3 for geno calling

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

# Check against shotgun geno
mismatchresults <- data.frame(Sample=popshotgun_selected$Samples,myBaits_match_shotgun=1,Twist_match_shotgun=2,het_geno=3,myBaits_match_shotgun_per=4,Twist_match_shotgun_per=5,different_allele_per=6)
for (i in 1:length(popshotgun_selected$Samples)) {
  sample_name <- popshotgun_selected$Samples[i] 
  mismatched_rows <- which(differences[,sample_name], arr.ind = TRUE)
  mismatched_row_names <- rownames(differences)[mismatched_rows]
  
  # Check against shotgun
  check_shotgun <- intersect(mismatched_row_names, rownames(shotgungeno))
  twist2 <- twistsnp_common[check_shotgun,] %>%
    select(sample_name)
  mybaits2 <- mybaitssnp_common[check_shotgun,] %>%
    select(sample_name)
  geno <- shotgungeno[check_shotgun,] %>%
    select(sample_name)
  twist2 <- tibble::rownames_to_column(twist2, var = "rowname")
  mybaits2 <- tibble::rownames_to_column(mybaits2, var = "rowname")
  geno <- tibble::rownames_to_column(geno, var = "rowname")
  geno2 <- geno %>% separate(sample_name, into =c("allele1","allele2"),sep=1) %>%
    mutate(het=(allele1!=allele2)*1)
  het <- sum(geno2$het, na.rm = TRUE)/sum(!is.na(geno2$het))*100
  
  merged_df <- mybaits2 %>%
    full_join(twist2, by = "rowname") %>%
    full_join(geno2, by="rowname")
  colnames(merged_df) <- c("rowname","myBaits","Twist","shotgun_allele1", "shotgun_allele2", "het")
  
  mybaits_match_shotgun <- sum(merged_df$myBaits == merged_df$shotgun_allele1 | merged_df$myBaits == merged_df$shotgun_allele2, na.rm=TRUE) 
  mismatchresults[i,2] <- mybaits_match_shotgun
  twist_match_shotgun <- sum(merged_df$Twist == merged_df$shotgun_allele1 | merged_df$Twist == merged_df$shotgun_allele2, na.rm=TRUE) 
  mismatchresults[i,3] <- twist_match_shotgun
  mismatchresults[i,4] <- het
  mismatchresults[i,5] <- mybaits_match_shotgun / merged_df%>% drop_na() %>% nrow() * 100
  twist_match_shotgun <- sum(merged_df$Twist == merged_df$shotgun_allele1 | merged_df$Twist == merged_df$shotgun_allele2, na.rm=TRUE) 
  mismatchresults[i,6] <- twist_match_shotgun / merged_df%>% drop_na() %>% nrow() * 100
  mismatchresults[i,7] <- sum(merged_df$shotgun_allele1 != merged_df$myBaits & merged_df$shotgun_allele2 != merged_df$myBaits & merged_df$shotgun_allele1 != merged_df$Twist & merged_df$shotgun_allele2 != merged_df$Twist, na.rm=TRUE) /  merged_df%>% drop_na() %>% nrow() *100
}

mismatchresults2 <- melt(as.data.table(mismatchresults), id.vars=c(1,2,3,4), variable.name="Method", value.name="Matches to shotgun (%)")
m1 <- ggplot(data=mismatchresults2,aes(x=Sample,y=`Matches to shotgun (%)`, fill=Method)) + 
  geom_bar(stat="identity",position="dodge",show.legend = TRUE) + 
  geom_text(aes(label = format(`Matches to shotgun (%)`,digits=1,scientific=FALSE)),position = position_dodge(width = 0.9), vjust = -0.3, size=3) +
  geom_text(aes(y = -max(`Matches to shotgun (%)`)*0.05,label = paste("het=",format(het_geno,digits=3,scientific=FALSE),"%",sep="")),size=3.5) +
  theme_bw() +
  scale_fill_manual(values = c("#F8766D", "#00BFC4", "orange"), name = "Method", labels = c("myBaits", "Twist","No match")) +
  theme(legend.position="bottom") +
  theme(axis.text.x = element_text(face="bold")) +
  theme(axis.title.x=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 

grid.arrange(p5,arrangeGrob(venn.plot, m1, ncol = 2))

##### pseudohaploid genotype for allele bias #####
## BRW001 ##
ac <- shotgungeno %>% select(BRW001) %>%
  filter(BRW001 == "AC" | BRW001 == "CA") 
ag <- shotgungeno %>% select(BRW001) %>%
  filter(BRW001 == "AG" | BRW001 == "GA") 
at <- shotgungeno %>% select(BRW001) %>%
  filter(BRW001 == "AT" | BRW001 == "TA")  
cg <- shotgungeno %>% select(BRW001) %>%
  filter(BRW001 == "CG" | BRW001 == "GC") 
ct <- shotgungeno %>% select(BRW001) %>%
  filter(BRW001 == "CT" | BRW001 == "TC") 
gt <- shotgungeno %>% select(BRW001) %>%
  filter(BRW001 == "GT" | BRW001 == "TG")  
het_brw001 <- rbind(ac,ag,at,cg,ct,gt)

#mybaits
genotype_counts <- list(
  AC = c(A = sum(mybaitssnp[rownames(mybaitssnp) %in% rownames(ac), ]$BRW001 == "A", na.rm = TRUE),
         C = sum(mybaitssnp[rownames(mybaitssnp) %in% rownames(ac), ]$BRW001 == "C", na.rm = TRUE)),
  AG = c(A = sum(mybaitssnp[rownames(mybaitssnp) %in% rownames(ag), ]$BRW001 == "A", na.rm = TRUE),
         G = sum(mybaitssnp[rownames(mybaitssnp) %in% rownames(ag), ]$BRW001 == "G", na.rm = TRUE)),
  AT = c(A = sum(mybaitssnp[rownames(mybaitssnp) %in% rownames(at), ]$BRW001 == "A", na.rm = TRUE),
         T = sum(mybaitssnp[rownames(mybaitssnp) %in% rownames(at), ]$BRW001 == "T", na.rm = TRUE)),
  CG = c(C = sum(mybaitssnp[rownames(mybaitssnp) %in% rownames(cg), ]$BRW001 == "C", na.rm = TRUE),
         G = sum(mybaitssnp[rownames(mybaitssnp) %in% rownames(cg), ]$BRW001 == "G", na.rm = TRUE)),
  CT = c(C = sum(mybaitssnp[rownames(mybaitssnp) %in% rownames(ct), ]$BRW001 == "C", na.rm = TRUE),
         T = sum(mybaitssnp[rownames(mybaitssnp) %in% rownames(ct), ]$BRW001 == "T", na.rm = TRUE)),
  GT = c(G = sum(mybaitssnp[rownames(mybaitssnp) %in% rownames(gt), ]$BRW001 == "G", na.rm = TRUE),
         T = sum(mybaitssnp[rownames(mybaitssnp) %in% rownames(gt), ]$BRW001 == "T", na.rm = TRUE)))
chi_square_results <- lapply(genotype_counts, function(counts) {
  chisq.test(x = counts, p = c(0.5, 0.5))  # Test if the counts for the two alleles are 50:50
})
p_values <- sapply(chi_square_results, function(x) x$p.value)
p_adjusted <- p.adjust(p_values, method = "bonferroni")
print(p_adjusted)
 
#twist
genotype_counts <- list(
  AC = c(A = sum(twistsnp[rownames(twistsnp) %in% rownames(ac), ]$BRW001 == "A", na.rm = TRUE),
         C = sum(twistsnp[rownames(twistsnp) %in% rownames(ac), ]$BRW001 == "C", na.rm = TRUE)),
  AG = c(A = sum(twistsnp[rownames(twistsnp) %in% rownames(ag), ]$BRW001 == "A", na.rm = TRUE),
         G = sum(twistsnp[rownames(twistsnp) %in% rownames(ag), ]$BRW001 == "G", na.rm = TRUE)),
  AT = c(A = sum(twistsnp[rownames(twistsnp) %in% rownames(at), ]$BRW001 == "A", na.rm = TRUE),
         T = sum(twistsnp[rownames(twistsnp) %in% rownames(at), ]$BRW001 == "T", na.rm = TRUE)),
  CG = c(C = sum(twistsnp[rownames(twistsnp) %in% rownames(cg), ]$BRW001 == "C", na.rm = TRUE),
         G = sum(twistsnp[rownames(twistsnp) %in% rownames(cg), ]$BRW001 == "G", na.rm = TRUE)),
  CT = c(C = sum(twistsnp[rownames(twistsnp) %in% rownames(ct), ]$BRW001 == "C", na.rm = TRUE),
         T = sum(twistsnp[rownames(twistsnp) %in% rownames(ct), ]$BRW001 == "T", na.rm = TRUE)),
  GT = c(G = sum(twistsnp[rownames(twistsnp) %in% rownames(gt), ]$BRW001 == "G", na.rm = TRUE),
         T = sum(twistsnp[rownames(twistsnp) %in% rownames(gt), ]$BRW001 == "T", na.rm = TRUE)))

chi_square_results <- lapply(genotype_counts, function(counts) {
  chisq.test(x = counts, p = c(0.5, 0.5))  # Test if the counts for the two alleles are 50:50
})
p_values <- sapply(chi_square_results, function(x) x$p.value)
p_adjusted <- p.adjust(p_values, method = "bonferroni")
print(genotype_counts)
print(p_adjusted)

## DVT014 ##
ac <- shotgungeno %>% select(DVT014) %>%
  filter(DVT014 == "AC" | DVT014 == "CA") 
ag <- shotgungeno %>% select(DVT014) %>%
  filter(DVT014 == "AG" | DVT014 == "GA") 
at <- shotgungeno %>% select(DVT014) %>%
  filter(DVT014 == "AT" | DVT014 == "TA")  
cg <- shotgungeno %>% select(DVT014) %>%
  filter(DVT014 == "CG" | DVT014 == "GC") 
ct <- shotgungeno %>% select(DVT014) %>%
  filter(DVT014 == "CT" | DVT014 == "TC") 
gt <- shotgungeno %>% select(DVT014) %>%
  filter(DVT014 == "GT" | DVT014 == "TG")  
het_dvt014 <- rbind(ac,ag,at,cg,ct,gt)

#mybaits
genotype_counts <- list(
  AC = c(A = sum(mybaitssnp[rownames(mybaitssnp) %in% rownames(ac), ]$DVT014 == "A", na.rm = TRUE),
         C = sum(mybaitssnp[rownames(mybaitssnp) %in% rownames(ac), ]$DVT014 == "C", na.rm = TRUE)),
  AG = c(A = sum(mybaitssnp[rownames(mybaitssnp) %in% rownames(ag), ]$DVT014 == "A", na.rm = TRUE),
         G = sum(mybaitssnp[rownames(mybaitssnp) %in% rownames(ag), ]$DVT014 == "G", na.rm = TRUE)),
  AT = c(A = sum(mybaitssnp[rownames(mybaitssnp) %in% rownames(at), ]$DVT014 == "A", na.rm = TRUE),
         T = sum(mybaitssnp[rownames(mybaitssnp) %in% rownames(at), ]$DVT014 == "T", na.rm = TRUE)),
  CG = c(C = sum(mybaitssnp[rownames(mybaitssnp) %in% rownames(cg), ]$DVT014 == "C", na.rm = TRUE),
         G = sum(mybaitssnp[rownames(mybaitssnp) %in% rownames(cg), ]$DVT014 == "G", na.rm = TRUE)),
  CT = c(C = sum(mybaitssnp[rownames(mybaitssnp) %in% rownames(ct), ]$DVT014 == "C", na.rm = TRUE),
         T = sum(mybaitssnp[rownames(mybaitssnp) %in% rownames(ct), ]$DVT014 == "T", na.rm = TRUE)),
  GT = c(G = sum(mybaitssnp[rownames(mybaitssnp) %in% rownames(gt), ]$DVT014 == "G", na.rm = TRUE),
         T = sum(mybaitssnp[rownames(mybaitssnp) %in% rownames(gt), ]$DVT014 == "T", na.rm = TRUE)))
chi_square_results <- lapply(genotype_counts, function(counts) {
  chisq.test(x = counts, p = c(0.5, 0.5))  # Test if the counts for the two alleles are 50:50
})
p_values <- sapply(chi_square_results, function(x) x$p.value)
p_adjusted <- p.adjust(p_values, method = "bonferroni")
print(genotype_counts)
print(p_adjusted)

#twist
genotype_counts <- list(
  AC = c(A = sum(twistsnp[rownames(twistsnp) %in% rownames(ac), ]$DVT014 == "A", na.rm = TRUE),
         C = sum(twistsnp[rownames(twistsnp) %in% rownames(ac), ]$DVT014 == "C", na.rm = TRUE)),
  AG = c(A = sum(twistsnp[rownames(twistsnp) %in% rownames(ag), ]$DVT014 == "A", na.rm = TRUE),
         G = sum(twistsnp[rownames(twistsnp) %in% rownames(ag), ]$DVT014 == "G", na.rm = TRUE)),
  AT = c(A = sum(twistsnp[rownames(twistsnp) %in% rownames(at), ]$DVT014 == "A", na.rm = TRUE),
         T = sum(twistsnp[rownames(twistsnp) %in% rownames(at), ]$DVT014 == "T", na.rm = TRUE)),
  CG = c(C = sum(twistsnp[rownames(twistsnp) %in% rownames(cg), ]$DVT014 == "C", na.rm = TRUE),
         G = sum(twistsnp[rownames(twistsnp) %in% rownames(cg), ]$DVT014 == "G", na.rm = TRUE)),
  CT = c(C = sum(twistsnp[rownames(twistsnp) %in% rownames(ct), ]$DVT014 == "C", na.rm = TRUE),
         T = sum(twistsnp[rownames(twistsnp) %in% rownames(ct), ]$DVT014 == "T", na.rm = TRUE)),
  GT = c(G = sum(twistsnp[rownames(twistsnp) %in% rownames(gt), ]$DVT014 == "G", na.rm = TRUE),
         T = sum(twistsnp[rownames(twistsnp) %in% rownames(gt), ]$DVT014 == "T", na.rm = TRUE)))

chi_square_results <- lapply(genotype_counts, function(counts) {
  chisq.test(x = counts, p = c(0.5, 0.5))  # Test if the counts for the two alleles are 50:50
})
p_values <- sapply(chi_square_results, function(x) x$p.value)
p_adjusted <- p.adjust(p_values, method = "bonferroni")
print(genotype_counts)
print(p_adjusted)

## NCP002 ##
ac <- shotgungeno %>% select(NCP002) %>%
  filter(NCP002 == "AC" | NCP002 == "CA") 
ag <- shotgungeno %>% select(NCP002) %>%
  filter(NCP002 == "AG" | NCP002 == "GA") 
at <- shotgungeno %>% select(NCP002) %>%
  filter(NCP002 == "AT" | NCP002 == "TA")  
cg <- shotgungeno %>% select(NCP002) %>%
  filter(NCP002 == "CG" | NCP002 == "GC") 
ct <- shotgungeno %>% select(NCP002) %>%
  filter(NCP002 == "CT" | NCP002 == "TC") 
gt <- shotgungeno %>% select(NCP002) %>%
  filter(NCP002 == "GT" | NCP002 == "TG")  
het_ncp002 <- rbind(ac,ag,at,cg,ct,gt)

#mybaits
genotype_counts <- list(
  AC = c(A = sum(mybaitssnp[rownames(mybaitssnp) %in% rownames(ac), ]$NCP002 == "A", na.rm = TRUE),
         C = sum(mybaitssnp[rownames(mybaitssnp) %in% rownames(ac), ]$NCP002 == "C", na.rm = TRUE)),
  AG = c(A = sum(mybaitssnp[rownames(mybaitssnp) %in% rownames(ag), ]$NCP002 == "A", na.rm = TRUE),
         G = sum(mybaitssnp[rownames(mybaitssnp) %in% rownames(ag), ]$NCP002 == "G", na.rm = TRUE)),
  AT = c(A = sum(mybaitssnp[rownames(mybaitssnp) %in% rownames(at), ]$NCP002 == "A", na.rm = TRUE),
         T = sum(mybaitssnp[rownames(mybaitssnp) %in% rownames(at), ]$NCP002 == "T", na.rm = TRUE)),
  CG = c(C = sum(mybaitssnp[rownames(mybaitssnp) %in% rownames(cg), ]$NCP002 == "C", na.rm = TRUE),
         G = sum(mybaitssnp[rownames(mybaitssnp) %in% rownames(cg), ]$NCP002 == "G", na.rm = TRUE)),
  CT = c(C = sum(mybaitssnp[rownames(mybaitssnp) %in% rownames(ct), ]$NCP002 == "C", na.rm = TRUE),
         T = sum(mybaitssnp[rownames(mybaitssnp) %in% rownames(ct), ]$NCP002 == "T", na.rm = TRUE)),
  GT = c(G = sum(mybaitssnp[rownames(mybaitssnp) %in% rownames(gt), ]$NCP002 == "G", na.rm = TRUE),
         T = sum(mybaitssnp[rownames(mybaitssnp) %in% rownames(gt), ]$NCP002 == "T", na.rm = TRUE)))
chi_square_results <- lapply(genotype_counts, function(counts) {
  chisq.test(x = counts, p = c(0.5, 0.5))  # Test if the counts for the two alleles are 50:50
})
p_values <- sapply(chi_square_results, function(x) x$p.value)
p_adjusted <- p.adjust(p_values, method = "bonferroni")
print(genotype_counts)
print(p_adjusted)

#twist
genotype_counts <- list(
  AC = c(A = sum(twistsnp[rownames(twistsnp) %in% rownames(ac), ]$NCP002 == "A", na.rm = TRUE),
         C = sum(twistsnp[rownames(twistsnp) %in% rownames(ac), ]$NCP002 == "C", na.rm = TRUE)),
  AG = c(A = sum(twistsnp[rownames(twistsnp) %in% rownames(ag), ]$NCP002 == "A", na.rm = TRUE),
         G = sum(twistsnp[rownames(twistsnp) %in% rownames(ag), ]$NCP002 == "G", na.rm = TRUE)),
  AT = c(A = sum(twistsnp[rownames(twistsnp) %in% rownames(at), ]$NCP002 == "A", na.rm = TRUE),
         T = sum(twistsnp[rownames(twistsnp) %in% rownames(at), ]$NCP002 == "T", na.rm = TRUE)),
  CG = c(C = sum(twistsnp[rownames(twistsnp) %in% rownames(cg), ]$NCP002 == "C", na.rm = TRUE),
         G = sum(twistsnp[rownames(twistsnp) %in% rownames(cg), ]$NCP002 == "G", na.rm = TRUE)),
  CT = c(C = sum(twistsnp[rownames(twistsnp) %in% rownames(ct), ]$NCP002 == "C", na.rm = TRUE),
         T = sum(twistsnp[rownames(twistsnp) %in% rownames(ct), ]$NCP002 == "T", na.rm = TRUE)),
  GT = c(G = sum(twistsnp[rownames(twistsnp) %in% rownames(gt), ]$NCP002 == "G", na.rm = TRUE),
         T = sum(twistsnp[rownames(twistsnp) %in% rownames(gt), ]$NCP002 == "T", na.rm = TRUE)))

chi_square_results <- lapply(genotype_counts, function(counts) {
  chisq.test(x = counts, p = c(0.5, 0.5))  # Test if the counts for the two alleles are 50:50
})
p_values <- sapply(chi_square_results, function(x) x$p.value)
p_adjusted <- p.adjust(p_values, method = "bonferroni")
print(genotype_counts)
print(p_adjusted)

## VKP001 ##
ac <- shotgungeno %>% select(VKP001) %>%
  filter(VKP001 == "AC" | VKP001 == "CA") 
ag <- shotgungeno %>% select(VKP001) %>%
  filter(VKP001 == "AG" | VKP001 == "GA") 
at <- shotgungeno %>% select(VKP001) %>%
  filter(VKP001 == "AT" | VKP001 == "TA")  
cg <- shotgungeno %>% select(VKP001) %>%
  filter(VKP001 == "CG" | VKP001 == "GC") 
ct <- shotgungeno %>% select(VKP001) %>%
  filter(VKP001 == "CT" | VKP001 == "TC") 
gt <- shotgungeno %>% select(VKP001) %>%
  filter(VKP001 == "GT" | VKP001 == "TG")  
het_vkp001 <- rbind(ac,ag,at,cg,ct,gt)

#mybaits
genotype_counts <- list(
  AC = c(A = sum(mybaitssnp[rownames(mybaitssnp) %in% rownames(ac), ]$VKP001 == "A", na.rm = TRUE),
         C = sum(mybaitssnp[rownames(mybaitssnp) %in% rownames(ac), ]$VKP001 == "C", na.rm = TRUE)),
  AG = c(A = sum(mybaitssnp[rownames(mybaitssnp) %in% rownames(ag), ]$VKP001 == "A", na.rm = TRUE),
         G = sum(mybaitssnp[rownames(mybaitssnp) %in% rownames(ag), ]$VKP001 == "G", na.rm = TRUE)),
  AT = c(A = sum(mybaitssnp[rownames(mybaitssnp) %in% rownames(at), ]$VKP001 == "A", na.rm = TRUE),
         T = sum(mybaitssnp[rownames(mybaitssnp) %in% rownames(at), ]$VKP001 == "T", na.rm = TRUE)),
  CG = c(C = sum(mybaitssnp[rownames(mybaitssnp) %in% rownames(cg), ]$VKP001 == "C", na.rm = TRUE),
         G = sum(mybaitssnp[rownames(mybaitssnp) %in% rownames(cg), ]$VKP001 == "G", na.rm = TRUE)),
  CT = c(C = sum(mybaitssnp[rownames(mybaitssnp) %in% rownames(ct), ]$VKP001 == "C", na.rm = TRUE),
         T = sum(mybaitssnp[rownames(mybaitssnp) %in% rownames(ct), ]$VKP001 == "T", na.rm = TRUE)),
  GT = c(G = sum(mybaitssnp[rownames(mybaitssnp) %in% rownames(gt), ]$VKP001 == "G", na.rm = TRUE),
         T = sum(mybaitssnp[rownames(mybaitssnp) %in% rownames(gt), ]$VKP001 == "T", na.rm = TRUE)))
chi_square_results <- lapply(genotype_counts, function(counts) {
  chisq.test(x = counts, p = c(0.5, 0.5))  # Test if the counts for the two alleles are 50:50
})
p_values <- sapply(chi_square_results, function(x) x$p.value)
p_adjusted <- p.adjust(p_values, method = "bonferroni")
print(genotype_counts)
print(p_adjusted)

#twist
genotype_counts <- list(
  AC = c(A = sum(twistsnp[rownames(twistsnp) %in% rownames(ac), ]$VKP001 == "A", na.rm = TRUE),
         C = sum(twistsnp[rownames(twistsnp) %in% rownames(ac), ]$VKP001 == "C", na.rm = TRUE)),
  AG = c(A = sum(twistsnp[rownames(twistsnp) %in% rownames(ag), ]$VKP001 == "A", na.rm = TRUE),
         G = sum(twistsnp[rownames(twistsnp) %in% rownames(ag), ]$VKP001 == "G", na.rm = TRUE)),
  AT = c(A = sum(twistsnp[rownames(twistsnp) %in% rownames(at), ]$VKP001 == "A", na.rm = TRUE),
         T = sum(twistsnp[rownames(twistsnp) %in% rownames(at), ]$VKP001 == "T", na.rm = TRUE)),
  CG = c(C = sum(twistsnp[rownames(twistsnp) %in% rownames(cg), ]$VKP001 == "C", na.rm = TRUE),
         G = sum(twistsnp[rownames(twistsnp) %in% rownames(cg), ]$VKP001 == "G", na.rm = TRUE)),
  CT = c(C = sum(twistsnp[rownames(twistsnp) %in% rownames(ct), ]$VKP001 == "C", na.rm = TRUE),
         T = sum(twistsnp[rownames(twistsnp) %in% rownames(ct), ]$VKP001 == "T", na.rm = TRUE)),
  GT = c(G = sum(twistsnp[rownames(twistsnp) %in% rownames(gt), ]$VKP001 == "G", na.rm = TRUE),
         T = sum(twistsnp[rownames(twistsnp) %in% rownames(gt), ]$VKP001 == "T", na.rm = TRUE)))

chi_square_results <- lapply(genotype_counts, function(counts) {
  chisq.test(x = counts, p = c(0.5, 0.5))  # Test if the counts for the two alleles are 50:50
})
p_values <- sapply(chi_square_results, function(x) x$p.value)
p_adjusted <- p.adjust(p_values, method = "bonferroni")
print(genotype_counts)
print(p_adjusted)

##### Allele depth for allele bias #####
# Table 1  and Supp figure 3 
#ind0:BRW001, ind1:DVT014, ind2:NCP002, ind3:VKP001
mybaitspos <- read.table(gzfile("mybaits_104kpanel_geno_maxmis_q20_dp3.pos.gz"), sep = "\t", header=TRUE)
twistpos <- read.table(gzfile("twist_104kpanel_geno_maxmis_q20_dp3.pos.gz"), sep = "\t", header=TRUE)
shotgunpos <- read.table(gzfile("shotgun_noudg_104kpanel_geno_maxmis_q20_dp3.pos.gz"), sep = "\t", header=TRUE)
mybaitscount <- read.table(gzfile("mybaits_104kpanel_geno_maxmis_q20_dp3.counts.gz"), sep = "\t", header=TRUE)
twistcount <- read.table(gzfile("twist_104kpanel_geno_maxmis_q20_dp3.counts.gz"), sep = "\t", header=TRUE)
shotguncount <- read.table(gzfile("shotgun_noudg_104kpanel_geno_maxmis_q20_dp3.counts.gz"), sep = "\t", header=TRUE)


mybaitspos <- read.table(gzfile("mybaits_104kpanel_geno_maxmis_q20_dp3.pos.gz"), sep = "\t", header=TRUE)
twistpos <- read.table(gzfile("twist_104kpanel_geno_maxmis_q20_dp3.pos.gz"), sep = "\t", header=TRUE)
shotgunpos <- read.table(gzfile("shotgun_noudg_104kpanel_geno_maxmis_q20_dp3.pos.gz"), sep = "\t", header=TRUE)
mybaitscount <- read.table(gzfile("mybaits_104kpanel_geno_maxmis_q20_dp3.counts.gz"), sep = "\t", header=TRUE)
twistcount <- read.table(gzfile("twist_104kpanel_geno_maxmis_q20_dp3.counts.gz"), sep = "\t", header=TRUE)
shotguncount <- read.table(gzfile("shotgun_noudg_104kpanel_geno_maxmis_q20_dp3.counts.gz"), sep = "\t", header=TRUE)

#Sample BRW001
mybaitsAD <- cbind(mybaitspos, mybaitscount[1:4]) %>%
  rename(A=4, C=5, G=6, T=7) %>%
  mutate(totDepth=A+C+G+T)  %>%
  mutate(scaff = paste(chr, pos, sep = "_")) %>%
  filter(totDepth >2) %>%
  rowwise() %>%
  filter(sum(c(A > 0, C > 0, G > 0, T > 0)) == 2) %>%
  mutate(
    AC = ifelse(A > 0 & C > 0, A/(A+C), 0),
    AG = ifelse(A > 0 & G > 0, A/(A+G), 0),
    AT = ifelse(A > 0 & T > 0, A/(A+T), 0),
    CG = ifelse(C > 0 & G > 0, C/(C+G), 0),
    CT = ifelse(C > 0 & T > 0, C/(C+T), 0),
    GT = ifelse(G > 0 & T > 0, G/(G+T), 0))
twistAD <- cbind(twistpos, twistcount[1:4]) %>%
  rename(A=4, C=5, G=6, T=7) %>%
  mutate(totDepth=A+C+G+T)  %>%
  mutate(scaff = paste(chr, pos, sep = "_")) %>%
  filter(totDepth >2) %>%
  rowwise() %>%
  filter(sum(c(A > 0, C > 0, G > 0, T > 0)) == 2) %>%
  mutate(
    AC = ifelse(A > 0 & C > 0, A/(A+C), 0),
    AG = ifelse(A > 0 & G > 0, A/(A+G), 0),
    AT = ifelse(A > 0 & T > 0, A/(A+T), 0),
    CG = ifelse(C > 0 & G > 0, C/(C+G), 0),
    CT = ifelse(C > 0 & T > 0, C/(C+T), 0),
    GT = ifelse(G > 0 & T > 0, G/(G+T), 0))
shotgunAD <- cbind(shotgunpos, shotguncount[1:4]) %>%
  rename(A=4, C=5, G=6, T=7) %>%
  mutate(totDepth=A+C+G+T)  %>%
  mutate(scaff = paste(chr, pos, sep = "_")) %>%
  filter(totDepth >2) %>%
  rowwise() %>%
  filter(sum(c(A > 0, C > 0, G > 0, T > 0)) == 2) %>%
  mutate(
    AC = ifelse(A > 0 & C > 0, A/(A+C), 0),
    AG = ifelse(A > 0 & G > 0, A/(A+G), 0),
    AT = ifelse(A > 0 & T > 0, A/(A+T), 0),
    CG = ifelse(C > 0 & G > 0, C/(C+G), 0),
    CT = ifelse(C > 0 & T > 0, C/(C+T), 0),
    GT = ifelse(G > 0 & T > 0, G/(G+T), 0))
common_rows <- Reduce(intersect, list(rownames(het_brw001), mybaitsAD$scaff, twistAD$scaff, shotgunAD$scaff)) #10407
mybaitsAD_filtered <- mybaitsAD[mybaitsAD$scaff %in% common_rows, ]
twistAD_filtered <- twistAD[twistAD$scaff %in% common_rows, ]
shotgunAD_filtered <- shotgunAD[shotgunAD$scaff %in% common_rows, ]

res_mybaits <- data.frame(Genotype=character(), estADratio=numeric(), p.value=numeric(), 
                          CI.low=numeric(), CI.high=numeric(), stringsAsFactors=FALSE)
logistic<-function(y=0){1/(1+exp(-y))}
for (i in c("AC", "AG", "AT", "CG", "CT", "GT")) {
  allele1 <- substr(i,1,1)
  allele2 <- substr(i,2,2)
  filtered_data <- mybaitsAD_filtered %>% filter(.data[[i]] > 0)
  allele1_counts <- filtered_data[[allele1]]
  allele2_counts <- filtered_data[[allele2]]
  lm1 <- glm(cbind(allele1_counts, allele2_counts)~1,data=filtered_data,family="quasibinomial")
  p.value<-summary(lm1)$coefficients[4]
  p.trans<-logistic(lm1$coefficients)
  CIlow.trans<-logistic(lm1$coefficients-summary(lm1)$coefficients[2]*1.96)
  CIup.trans<-logistic(lm1$coefficients+summary(lm1)$coefficients[2]*1.96)
  res_mybaits <- rbind(res_mybaits, data.frame(Genotype=i, estADratio=p.trans, p.value=p.value, 
                                               CI.low=CIlow.trans, CI.high=CIup.trans)) }

res_twist <- data.frame(Genotype=character(), estADratio=numeric(), p.value=numeric(), 
                        CI.low=numeric(), CI.high=numeric(), stringsAsFactors=FALSE)
logistic<-function(y=0){1/(1+exp(-y))}
for (i in c("AC", "AG", "AT", "CG", "CT", "GT")) {
  allele1 <- substr(i,1,1)
  allele2 <- substr(i,2,2)
  filtered_data <- twistAD_filtered %>% filter(.data[[i]] > 0)
  allele1_counts <- filtered_data[[allele1]]
  allele2_counts <- filtered_data[[allele2]]
  lm1 <- glm(cbind(allele1_counts, allele2_counts)~1,data=filtered_data,family="quasibinomial")
  p.value<-summary(lm1)$coefficients[4]
  p.trans<-logistic(lm1$coefficients)
  CIlow.trans<-logistic(lm1$coefficients-summary(lm1)$coefficients[2]*1.96)
  CIup.trans<-logistic(lm1$coefficients+summary(lm1)$coefficients[2]*1.96)
  res_twist <- rbind(res_twist, data.frame(Genotype=i, estADratio=p.trans, p.value=p.value, 
                                           CI.low=CIlow.trans, CI.high=CIup.trans)) }

res_shotgun <- data.frame(Genotype=character(), estADratio=numeric(), p.value=numeric(), 
                          CI.low=numeric(), CI.high=numeric(), stringsAsFactors=FALSE)
logistic<-function(y=0){1/(1+exp(-y))}
for (i in c("AC", "AG", "AT", "CG", "CT", "GT")) {
  allele1 <- substr(i,1,1)
  allele2 <- substr(i,2,2)
  filtered_data <- shotgunAD_filtered %>% filter(.data[[i]] > 0)
  allele1_counts <- filtered_data[[allele1]]
  allele2_counts <- filtered_data[[allele2]]
  lm1 <- glm(cbind(allele1_counts, allele2_counts)~1,data=filtered_data,family="quasibinomial")
  p.value<-summary(lm1)$coefficients[4]
  p.trans<-logistic(lm1$coefficients)
  CIlow.trans<-logistic(lm1$coefficients-summary(lm1)$coefficients[2]*1.96)
  CIup.trans<-logistic(lm1$coefficients+summary(lm1)$coefficients[2]*1.96)
  res_shotgun <- rbind(res_shotgun, data.frame(Genotype=i, estADratio=p.trans, p.value=p.value, 
                                               CI.low=CIlow.trans, CI.high=CIup.trans)) }
res_shotgun
res_mybaits
res_twist

#### VIOLIN PLOT OF ALLELE BIAS
df1 <- data.frame(Frequency = shotgunAD_filtered%>%filter(AC>0)%>%select(AC), Group = "shotgun")
df2 <- data.frame(Frequency = mybaitsAD_filtered%>%filter(AC>0)%>%select(AC), Group = "myBaits")
df3 <- data.frame(Frequency = twistAD_filtered%>%filter(AC>0)%>%select(AC), Group = "Twist")
combined_AC <- rbind(df1, df2, df3) %>% rename("Frequency"=1) %>%
  mutate(Group = factor(Group, levels = c("shotgun", "myBaits", "Twist"))) 
df1 <- data.frame(Frequency = shotgunAD_filtered%>%filter(AG>0)%>%select(AG), Group = "shotgun")
df2 <- data.frame(Frequency = mybaitsAD_filtered%>%filter(AG>0)%>%select(AG), Group = "myBaits")
df3 <- data.frame(Frequency = twistAD_filtered%>%filter(AG>0)%>%select(AG), Group = "Twist")
combined_AG <- rbind(df1, df2, df3) %>% rename("Frequency"=1) %>%
  mutate(Group = factor(Group, levels = c("shotgun", "myBaits", "Twist"))) 
df1 <- data.frame(Frequency = shotgunAD_filtered%>%filter(AT>0)%>%select(AT), Group = "shotgun")
df2 <- data.frame(Frequency = mybaitsAD_filtered%>%filter(AT>0)%>%select(AT), Group = "myBaits")
df3 <- data.frame(Frequency = twistAD_filtered%>%filter(AT>0)%>%select(AT), Group = "Twist")
combined_AT <- rbind(df1, df2, df3) %>% rename("Frequency"=1) %>%
  mutate(Group = factor(Group, levels = c("shotgun", "myBaits", "Twist"))) 
df1 <- data.frame(Frequency = shotgunAD_filtered%>%filter(CG>0)%>%select(CG), Group = "shotgun")
df2 <- data.frame(Frequency = mybaitsAD_filtered%>%filter(CG>0)%>%select(CG), Group = "myBaits")
df3 <- data.frame(Frequency = twistAD_filtered%>%filter(CG>0)%>%select(CG), Group = "Twist")
combined_CG <- rbind(df1, df2, df3) %>% rename("Frequency"=1) %>%
  mutate(Group = factor(Group, levels = c("shotgun", "myBaits", "Twist"))) 
df1 <- data.frame(Frequency = shotgunAD_filtered%>%filter(CT>0)%>%select(CT), Group = "shotgun")
df2 <- data.frame(Frequency = mybaitsAD_filtered%>%filter(CT>0)%>%select(CT), Group = "myBaits")
df3 <- data.frame(Frequency = twistAD_filtered%>%filter(CT>0)%>%select(CT), Group = "Twist")
combined_CT <- rbind(df1, df2, df3) %>% rename("Frequency"=1) %>%
  mutate(Group = factor(Group, levels = c("shotgun", "myBaits", "Twist"))) 
df1 <- data.frame(Frequency = shotgunAD_filtered%>%filter(GT>0)%>%select(GT), Group = "shotgun")
df2 <- data.frame(Frequency = mybaitsAD_filtered%>%filter(GT>0)%>%select(GT), Group = "myBaits")
df3 <- data.frame(Frequency = twistAD_filtered%>%filter(GT>0)%>%select(GT), Group = "Twist")
combined_GT <- rbind(df1, df2, df3) %>% rename("Frequency"=1) %>%
  mutate(Group = factor(Group, levels = c("shotgun", "myBaits", "Twist"))) 
a1 <- ggplot(combined_AG, aes(x = Group, y = Frequency, fill = Group)) +
  geom_violin(trim = TRUE) +  # Use trim = TRUE if you want to trim the tails
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +  # Add boxplot inside violin
  scale_fill_manual(values = c("orange","#F8766D", "#00BFC4"), name = "Method", labels = c("shotgun","myBaits", "Twist")) +
  labs(y="BRW001") +
  theme_minimal() +
  theme(legend.position = "none") + 
  theme(axis.title.x = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),axis.line.x = element_blank())
a2 <- ggplot(combined_CT, aes(x = Group, y = Frequency, fill = Group)) +
  geom_violin(trim = TRUE) +  # Use trim = TRUE if you want to trim the tails
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +  
  scale_fill_manual(values = c("orange","#F8766D", "#00BFC4"), name = "Method", labels = c("shotgun","myBaits", "Twist")) +
  theme_minimal() +
  theme(legend.position = "none") + 
  theme(axis.title.x = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),axis.line.x = element_blank(),axis.title.y = element_blank())
a3 <- ggplot(combined_AC, aes(x = Group, y = Frequency, fill = Group)) +
  geom_violin(trim = TRUE) +  # Use trim = TRUE if you want to trim the tails
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +  
  scale_fill_manual(values = c("orange","#F8766D", "#00BFC4"), name = "Method", labels = c("shotgun","myBaits", "Twist")) +
  theme_minimal() +
  theme(legend.position = "none") + 
  theme(axis.title.x = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),axis.line.x = element_blank(),axis.title.y = element_blank())
a4 <- ggplot(combined_AT, aes(x = Group, y = Frequency, fill = Group)) +
  geom_violin(trim = TRUE) +  # Use trim = TRUE if you want to trim the tails
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +  
  scale_fill_manual(values = c("orange","#F8766D", "#00BFC4"), name = "Method", labels = c("shotgun","myBaits", "Twist")) +
  theme_minimal() +
  theme(legend.position = "none") + 
  theme(axis.title.x = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),axis.line.x = element_blank(),axis.title.y = element_blank())
a5 <- ggplot(combined_CG, aes(x = Group, y = Frequency, fill = Group)) +
  geom_violin(trim = TRUE) +  # Use trim = TRUE if you want to trim the tails
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +  
  scale_fill_manual(values = c("orange","#F8766D", "#00BFC4"), name = "Method", labels = c("shotgun","myBaits", "Twist")) +
  theme_minimal() +
  theme(legend.position = "none") + 
  theme(axis.title.x = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),axis.line.x = element_blank(),axis.title.y = element_blank())
a6 <- ggplot(combined_GT, aes(x = Group, y = Frequency, fill = Group)) +
  geom_violin(trim = TRUE) +  # Use trim = TRUE if you want to trim the tails
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +  
  scale_fill_manual(values = c("orange","#F8766D", "#00BFC4"), name = "Method", labels = c("Shotgun","myBaits", "Twist")) +
  theme_minimal() +
  theme(legend.position = "none") + 
  theme(axis.title.x = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),axis.line.x = element_blank(),axis.title.y = element_blank())
# Create a function to extract the legend
get_legend <- function(my_plot) {
  g <- ggplotGrob(my_plot)
  legend <- g$grobs[[which(sapply(g$grobs, function(x) x$name) == "guide-box")]]
  if (is.null(legend)) stop("Legend not found in the plot!")
  return(legend)}
shared_legend <- get_legend(
  ggplot(combined_GT, aes(x = Group, y = Frequency, fill = Group)) +
    geom_violin(trim = TRUE) +  # Use trim = TRUE if you want to trim the tails
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +  
    scale_fill_manual(values = c("orange","#F8766D", "#00BFC4"), name = "Method", labels = c("Shotgun","myBaits", "Twist")) +
    theme_minimal() +
    theme(legend.position = "top") )

#Sample DVT014
mybaitsAD <- cbind(mybaitspos, mybaitscount[5:8]) %>%
  rename(A=4, C=5, G=6, T=7) %>%
  mutate(totDepth=A+C+G+T)  %>%
  mutate(scaff = paste(chr, pos, sep = "_")) %>%
  filter(totDepth >2) %>%
  rowwise() %>%
  filter(sum(c(A > 0, C > 0, G > 0, T > 0)) == 2) %>%
  mutate(
    AC = ifelse(A > 0 & C > 0, A/(A+C), 0),
    AG = ifelse(A > 0 & G > 0, A/(A+G), 0),
    AT = ifelse(A > 0 & T > 0, A/(A+T), 0),
    CG = ifelse(C > 0 & G > 0, C/(C+G), 0),
    CT = ifelse(C > 0 & T > 0, C/(C+T), 0),
    GT = ifelse(G > 0 & T > 0, G/(G+T), 0))
twistAD <- cbind(twistpos, twistcount[5:8]) %>%
  rename(A=4, C=5, G=6, T=7) %>%
  mutate(totDepth=A+C+G+T)  %>%
  mutate(scaff = paste(chr, pos, sep = "_")) %>%
  filter(totDepth >2) %>%
  rowwise() %>%
  filter(sum(c(A > 0, C > 0, G > 0, T > 0)) == 2) %>%
  mutate(
    AC = ifelse(A > 0 & C > 0, A/(A+C), 0),
    AG = ifelse(A > 0 & G > 0, A/(A+G), 0),
    AT = ifelse(A > 0 & T > 0, A/(A+T), 0),
    CG = ifelse(C > 0 & G > 0, C/(C+G), 0),
    CT = ifelse(C > 0 & T > 0, C/(C+T), 0),
    GT = ifelse(G > 0 & T > 0, G/(G+T), 0))
shotgunAD <- cbind(shotgunpos, shotguncount[5:8]) %>%
  rename(A=4, C=5, G=6, T=7) %>%
  mutate(totDepth=A+C+G+T)  %>%
  mutate(scaff = paste(chr, pos, sep = "_")) %>%
  filter(totDepth >2) %>%
  rowwise() %>%
  filter(sum(c(A > 0, C > 0, G > 0, T > 0)) == 2) %>%
  mutate(
    AC = ifelse(A > 0 & C > 0, A/(A+C), 0),
    AG = ifelse(A > 0 & G > 0, A/(A+G), 0),
    AT = ifelse(A > 0 & T > 0, A/(A+T), 0),
    CG = ifelse(C > 0 & G > 0, C/(C+G), 0),
    CT = ifelse(C > 0 & T > 0, C/(C+T), 0),
    GT = ifelse(G > 0 & T > 0, G/(G+T), 0))
common_rows <- Reduce(intersect, list(rownames(het_dvt014), mybaitsAD$scaff, twistAD$scaff, shotgunAD$scaff)) #3019
mybaitsAD_filtered <- mybaitsAD[mybaitsAD$scaff %in% common_rows, ]
twistAD_filtered <- twistAD[twistAD$scaff %in% common_rows, ]
shotgunAD_filtered <- shotgunAD[shotgunAD$scaff %in% common_rows, ]

res_mybaits <- data.frame(Genotype=character(), estADratio=numeric(), p.value=numeric(), 
                          CI.low=numeric(), CI.high=numeric(), stringsAsFactors=FALSE)
logistic<-function(y=0){1/(1+exp(-y))}
for (i in c("AC", "AG", "AT", "CG", "CT", "GT")) {
  allele1 <- substr(i,1,1)
  allele2 <- substr(i,2,2)
  filtered_data <- mybaitsAD_filtered %>% filter(.data[[i]] > 0)
  allele1_counts <- filtered_data[[allele1]]
  allele2_counts <- filtered_data[[allele2]]
  lm1 <- glm(cbind(allele1_counts, allele2_counts)~1,data=filtered_data,family="quasibinomial")
  p.value<-summary(lm1)$coefficients[4]
  p.trans<-logistic(lm1$coefficients)
  CIlow.trans<-logistic(lm1$coefficients-summary(lm1)$coefficients[2]*1.96)
  CIup.trans<-logistic(lm1$coefficients+summary(lm1)$coefficients[2]*1.96)
  res_mybaits <- rbind(res_mybaits, data.frame(Genotype=i, estADratio=p.trans, p.value=p.value, 
                                               CI.low=CIlow.trans, CI.high=CIup.trans)) }

res_twist <- data.frame(Genotype=character(), estADratio=numeric(), p.value=numeric(), 
                        CI.low=numeric(), CI.high=numeric(), stringsAsFactors=FALSE)
logistic<-function(y=0){1/(1+exp(-y))}
for (i in c("AC", "AG", "AT", "CG", "CT", "GT")) {
  allele1 <- substr(i,1,1)
  allele2 <- substr(i,2,2)
  filtered_data <- twistAD_filtered %>% filter(.data[[i]] > 0)
  allele1_counts <- filtered_data[[allele1]]
  allele2_counts <- filtered_data[[allele2]]
  lm1 <- glm(cbind(allele1_counts, allele2_counts)~1,data=filtered_data,family="quasibinomial")
  p.value<-summary(lm1)$coefficients[4]
  p.trans<-logistic(lm1$coefficients)
  CIlow.trans<-logistic(lm1$coefficients-summary(lm1)$coefficients[2]*1.96)
  CIup.trans<-logistic(lm1$coefficients+summary(lm1)$coefficients[2]*1.96)
  res_twist <- rbind(res_twist, data.frame(Genotype=i, estADratio=p.trans, p.value=p.value, 
                                           CI.low=CIlow.trans, CI.high=CIup.trans)) }

res_shotgun <- data.frame(Genotype=character(), estADratio=numeric(), p.value=numeric(), 
                          CI.low=numeric(), CI.high=numeric(), stringsAsFactors=FALSE)
logistic<-function(y=0){1/(1+exp(-y))}
for (i in c("AC", "AG", "AT", "CG", "CT", "GT")) {
  allele1 <- substr(i,1,1)
  allele2 <- substr(i,2,2)
  filtered_data <- shotgunAD_filtered %>% filter(.data[[i]] > 0)
  allele1_counts <- filtered_data[[allele1]]
  allele2_counts <- filtered_data[[allele2]]
  lm1 <- glm(cbind(allele1_counts, allele2_counts)~1,data=filtered_data,family="quasibinomial")
  p.value<-summary(lm1)$coefficients[4]
  p.trans<-logistic(lm1$coefficients)
  CIlow.trans<-logistic(lm1$coefficients-summary(lm1)$coefficients[2]*1.96)
  CIup.trans<-logistic(lm1$coefficients+summary(lm1)$coefficients[2]*1.96)
  res_shotgun <- rbind(res_shotgun, data.frame(Genotype=i, estADratio=p.trans, p.value=p.value, 
                                               CI.low=CIlow.trans, CI.high=CIup.trans)) }
res_shotgun
res_mybaits
res_twist

#### VIOLIN PLOT OF ALLELE BIAS
df1 <- data.frame(Frequency = shotgunAD_filtered%>%filter(AC>0)%>%select(AC), Group = "shotgun")
df2 <- data.frame(Frequency = mybaitsAD_filtered%>%filter(AC>0)%>%select(AC), Group = "myBaits")
df3 <- data.frame(Frequency = twistAD_filtered%>%filter(AC>0)%>%select(AC), Group = "Twist")
combined_AC <- rbind(df1, df2, df3) %>% rename("Frequency"=1) %>%
  mutate(Group = factor(Group, levels = c("shotgun", "myBaits", "Twist"))) 
df1 <- data.frame(Frequency = shotgunAD_filtered%>%filter(AG>0)%>%select(AG), Group = "shotgun")
df2 <- data.frame(Frequency = mybaitsAD_filtered%>%filter(AG>0)%>%select(AG), Group = "myBaits")
df3 <- data.frame(Frequency = twistAD_filtered%>%filter(AG>0)%>%select(AG), Group = "Twist")
combined_AG <- rbind(df1, df2, df3) %>% rename("Frequency"=1) %>%
  mutate(Group = factor(Group, levels = c("shotgun", "myBaits", "Twist"))) 
df1 <- data.frame(Frequency = shotgunAD_filtered%>%filter(AT>0)%>%select(AT), Group = "shotgun")
df2 <- data.frame(Frequency = mybaitsAD_filtered%>%filter(AT>0)%>%select(AT), Group = "myBaits")
df3 <- data.frame(Frequency = twistAD_filtered%>%filter(AT>0)%>%select(AT), Group = "Twist")
combined_AT <- rbind(df1, df2, df3) %>% rename("Frequency"=1) %>%
  mutate(Group = factor(Group, levels = c("shotgun", "myBaits", "Twist"))) 
df1 <- data.frame(Frequency = shotgunAD_filtered%>%filter(CG>0)%>%select(CG), Group = "shotgun")
df2 <- data.frame(Frequency = mybaitsAD_filtered%>%filter(CG>0)%>%select(CG), Group = "myBaits")
df3 <- data.frame(Frequency = twistAD_filtered%>%filter(CG>0)%>%select(CG), Group = "Twist")
combined_CG <- rbind(df1, df2, df3) %>% rename("Frequency"=1) %>%
  mutate(Group = factor(Group, levels = c("shotgun", "myBaits", "Twist"))) 
df1 <- data.frame(Frequency = shotgunAD_filtered%>%filter(CT>0)%>%select(CT), Group = "shotgun")
df2 <- data.frame(Frequency = mybaitsAD_filtered%>%filter(CT>0)%>%select(CT), Group = "myBaits")
df3 <- data.frame(Frequency = twistAD_filtered%>%filter(CT>0)%>%select(CT), Group = "Twist")
combined_CT <- rbind(df1, df2, df3) %>% rename("Frequency"=1) %>%
  mutate(Group = factor(Group, levels = c("shotgun", "myBaits", "Twist"))) 
df1 <- data.frame(Frequency = shotgunAD_filtered%>%filter(GT>0)%>%select(GT), Group = "shotgun")
df2 <- data.frame(Frequency = mybaitsAD_filtered%>%filter(GT>0)%>%select(GT), Group = "myBaits")
df3 <- data.frame(Frequency = twistAD_filtered%>%filter(GT>0)%>%select(GT), Group = "Twist")
combined_GT <- rbind(df1, df2, df3) %>% rename("Frequency"=1) %>%
  mutate(Group = factor(Group, levels = c("shotgun", "myBaits", "Twist"))) 
b1 <- ggplot(combined_AG, aes(x = Group, y = Frequency, fill = Group)) +
  geom_violin(trim = TRUE) +  # Use trim = TRUE if you want to trim the tails
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +  # Add boxplot inside violin
  scale_fill_manual(values = c("orange","#F8766D", "#00BFC4"), name = "Method", labels = c("shotgun","myBaits", "Twist")) +
  labs(y="DVT014") +
  theme_minimal() +
  theme(legend.position = "none") + 
  theme(axis.title.x = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),axis.line.x = element_blank())
b2 <- ggplot(combined_CT, aes(x = Group, y = Frequency, fill = Group)) +
  geom_violin(trim = TRUE) +  # Use trim = TRUE if you want to trim the tails
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +  
  scale_fill_manual(values = c("orange","#F8766D", "#00BFC4"), name = "Method", labels = c("shotgun","myBaits", "Twist")) +
  theme_minimal() +
  theme(legend.position = "none") + 
  theme(axis.title.x = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),axis.line.x = element_blank(),axis.title.y = element_blank())
b3 <- ggplot(combined_AC, aes(x = Group, y = Frequency, fill = Group)) +
  geom_violin(trim = TRUE) +  # Use trim = TRUE if you want to trim the tails
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +  
  scale_fill_manual(values = c("orange","#F8766D", "#00BFC4"), name = "Method", labels = c("shotgun","myBaits", "Twist")) +
  theme_minimal() +
  theme(legend.position = "none") + 
  theme(axis.title.x = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),axis.line.x = element_blank(),axis.title.y = element_blank())
b4 <- ggplot(combined_AT, aes(x = Group, y = Frequency, fill = Group)) +
  geom_violin(trim = TRUE) +  # Use trim = TRUE if you want to trim the tails
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +  
  scale_fill_manual(values = c("orange","#F8766D", "#00BFC4"), name = "Method", labels = c("shotgun","myBaits", "Twist")) +
  theme_minimal() +
  theme(legend.position = "none") + 
  theme(axis.title.x = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),axis.line.x = element_blank(),axis.title.y = element_blank())
b5 <- ggplot(combined_CG, aes(x = Group, y = Frequency, fill = Group)) +
  geom_violin(trim = TRUE) +  # Use trim = TRUE if you want to trim the tails
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +  
  scale_fill_manual(values = c("orange","#F8766D", "#00BFC4"), name = "Method", labels = c("shotgun","myBaits", "Twist")) +
  theme_minimal() +
  theme(legend.position = "none") + 
  theme(axis.title.x = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),axis.line.x = element_blank(),axis.title.y = element_blank())
b6 <- ggplot(combined_GT, aes(x = Group, y = Frequency, fill = Group)) +
  geom_violin(trim = TRUE) +  # Use trim = TRUE if you want to trim the tails
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +  
  scale_fill_manual(values = c("orange","#F8766D", "#00BFC4"), name = "Method", labels = c("shotgun","myBaits", "Twist")) +
  theme_minimal() +
  theme(legend.position = "none") + 
  theme(axis.title.x = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),axis.line.x = element_blank(),axis.title.y = element_blank())

#Sample NCP002
##mybaits
mybaitsAD <- cbind(mybaitspos, mybaitscount[9:12]) %>%
  rename(A=4, C=5, G=6, T=7) %>%
  mutate(totDepth=A+C+G+T)  %>%
  mutate(scaff = paste(chr, pos, sep = "_")) %>%
  filter(totDepth >2) %>%
  rowwise() %>%
  filter(sum(c(A > 0, C > 0, G > 0, T > 0)) == 2) %>%
  mutate(
    AC = ifelse(A > 0 & C > 0, A/(A+C), 0),
    AG = ifelse(A > 0 & G > 0, A/(A+G), 0),
    AT = ifelse(A > 0 & T > 0, A/(A+T), 0),
    CG = ifelse(C > 0 & G > 0, C/(C+G), 0),
    CT = ifelse(C > 0 & T > 0, C/(C+T), 0),
    GT = ifelse(G > 0 & T > 0, G/(G+T), 0))
twistAD <- cbind(twistpos, twistcount[9:12]) %>%
  rename(A=4, C=5, G=6, T=7) %>%
  mutate(totDepth=A+C+G+T)  %>%
  mutate(scaff = paste(chr, pos, sep = "_")) %>%
  filter(totDepth >2) %>%
  rowwise() %>%
  filter(sum(c(A > 0, C > 0, G > 0, T > 0)) == 2) %>%
  mutate(
    AC = ifelse(A > 0 & C > 0, A/(A+C), 0),
    AG = ifelse(A > 0 & G > 0, A/(A+G), 0),
    AT = ifelse(A > 0 & T > 0, A/(A+T), 0),
    CG = ifelse(C > 0 & G > 0, C/(C+G), 0),
    CT = ifelse(C > 0 & T > 0, C/(C+T), 0),
    GT = ifelse(G > 0 & T > 0, G/(G+T), 0))
shotgunAD <- cbind(shotgunpos, shotguncount[9:12]) %>%
  rename(A=4, C=5, G=6, T=7) %>%
  mutate(totDepth=A+C+G+T)  %>%
  mutate(scaff = paste(chr, pos, sep = "_")) %>%
  filter(totDepth >2) %>%
  rowwise() %>%
  filter(sum(c(A > 0, C > 0, G > 0, T > 0)) == 2) %>%
  mutate(
    AC = ifelse(A > 0 & C > 0, A/(A+C), 0),
    AG = ifelse(A > 0 & G > 0, A/(A+G), 0),
    AT = ifelse(A > 0 & T > 0, A/(A+T), 0),
    CG = ifelse(C > 0 & G > 0, C/(C+G), 0),
    CT = ifelse(C > 0 & T > 0, C/(C+T), 0),
    GT = ifelse(G > 0 & T > 0, G/(G+T), 0))
common_rows <- Reduce(intersect, list(rownames(het_ncp002), mybaitsAD$scaff, twistAD$scaff, shotgunAD$scaff)) #67
mybaitsAD_filtered <- mybaitsAD[mybaitsAD$scaff %in% common_rows, ]
twistAD_filtered <- twistAD[twistAD$scaff %in% common_rows, ]
shotgunAD_filtered <- shotgunAD[shotgunAD$scaff %in% common_rows, ]

res_mybaits <- data.frame(Genotype=character(), estADratio=numeric(), p.value=numeric(), 
                          CI.low=numeric(), CI.high=numeric(), stringsAsFactors=FALSE)
logistic<-function(y=0){1/(1+exp(-y))}
for (i in c("AC", "AG", "AT", "CG", "CT", "GT")) {
  allele1 <- substr(i,1,1)
  allele2 <- substr(i,2,2)
  filtered_data <- mybaitsAD_filtered %>% filter(.data[[i]] > 0)
  allele1_counts <- filtered_data[[allele1]]
  allele2_counts <- filtered_data[[allele2]]
  lm1 <- glm(cbind(allele1_counts, allele2_counts)~1,data=filtered_data,family="quasibinomial")
  p.value<-summary(lm1)$coefficients[4]
  p.trans<-logistic(lm1$coefficients)
  CIlow.trans<-logistic(lm1$coefficients-summary(lm1)$coefficients[2]*1.96)
  CIup.trans<-logistic(lm1$coefficients+summary(lm1)$coefficients[2]*1.96)
  res_mybaits <- rbind(res_mybaits, data.frame(Genotype=i, estADratio=p.trans, p.value=p.value, 
                                               CI.low=CIlow.trans, CI.high=CIup.trans)) }

res_twist <- data.frame(Genotype=character(), estADratio=numeric(), p.value=numeric(), 
                        CI.low=numeric(), CI.high=numeric(), stringsAsFactors=FALSE)
logistic<-function(y=0){1/(1+exp(-y))}
for (i in c("AC", "AG", "AT", "CG", "CT", "GT")) {
  allele1 <- substr(i,1,1)
  allele2 <- substr(i,2,2)
  filtered_data <- twistAD_filtered %>% filter(.data[[i]] > 0)
  allele1_counts <- filtered_data[[allele1]]
  allele2_counts <- filtered_data[[allele2]]
  lm1 <- glm(cbind(allele1_counts, allele2_counts)~1,data=filtered_data,family="quasibinomial")
  p.value<-summary(lm1)$coefficients[4]
  p.trans<-logistic(lm1$coefficients)
  CIlow.trans<-logistic(lm1$coefficients-summary(lm1)$coefficients[2]*1.96)
  CIup.trans<-logistic(lm1$coefficients+summary(lm1)$coefficients[2]*1.96)
  res_twist <- rbind(res_twist, data.frame(Genotype=i, estADratio=p.trans, p.value=p.value, 
                                           CI.low=CIlow.trans, CI.high=CIup.trans)) }

res_shotgun <- data.frame(Genotype=character(), estADratio=numeric(), p.value=numeric(), 
                          CI.low=numeric(), CI.high=numeric(), stringsAsFactors=FALSE)
logistic<-function(y=0){1/(1+exp(-y))}
for (i in c("AC", "AG", "AT", "CG", "CT", "GT")) {
  allele1 <- substr(i,1,1)
  allele2 <- substr(i,2,2)
  filtered_data <- shotgunAD_filtered %>% filter(.data[[i]] > 0)
  allele1_counts <- filtered_data[[allele1]]
  allele2_counts <- filtered_data[[allele2]]
  lm1 <- glm(cbind(allele1_counts, allele2_counts)~1,data=filtered_data,family="quasibinomial")
  p.value<-summary(lm1)$coefficients[4]
  p.trans<-logistic(lm1$coefficients)
  CIlow.trans<-logistic(lm1$coefficients-summary(lm1)$coefficients[2]*1.96)
  CIup.trans<-logistic(lm1$coefficients+summary(lm1)$coefficients[2]*1.96)
  res_shotgun <- rbind(res_shotgun, data.frame(Genotype=i, estADratio=p.trans, p.value=p.value, 
                                               CI.low=CIlow.trans, CI.high=CIup.trans)) }
res_shotgun
res_mybaits
res_twist

#### VIOLIN PLOT OF ALLELE BIAS
df1 <- data.frame(Frequency = shotgunAD_filtered%>%filter(AC>0)%>%select(AC), Group = "shotgun")
df2 <- data.frame(Frequency = mybaitsAD_filtered%>%filter(AC>0)%>%select(AC), Group = "myBaits")
df3 <- data.frame(Frequency = twistAD_filtered%>%filter(AC>0)%>%select(AC), Group = "Twist")
combined_AC <- rbind(df1, df2, df3) %>% rename("Frequency"=1) %>%
  mutate(Group = factor(Group, levels = c("shotgun", "myBaits", "Twist"))) 
df1 <- data.frame(Frequency = shotgunAD_filtered%>%filter(AG>0)%>%select(AG), Group = "shotgun")
df2 <- data.frame(Frequency = mybaitsAD_filtered%>%filter(AG>0)%>%select(AG), Group = "myBaits")
df3 <- data.frame(Frequency = twistAD_filtered%>%filter(AG>0)%>%select(AG), Group = "Twist")
combined_AG <- rbind(df1, df2, df3) %>% rename("Frequency"=1) %>%
  mutate(Group = factor(Group, levels = c("shotgun", "myBaits", "Twist"))) 
df1 <- data.frame(Frequency = shotgunAD_filtered%>%filter(AT>0)%>%select(AT), Group = "shotgun")
df2 <- data.frame(Frequency = mybaitsAD_filtered%>%filter(AT>0)%>%select(AT), Group = "myBaits")
df3 <- data.frame(Frequency = twistAD_filtered%>%filter(AT>0)%>%select(AT), Group = "Twist")
combined_AT <- rbind(df1, df2, df3) %>% rename("Frequency"=1) %>%
  mutate(Group = factor(Group, levels = c("shotgun", "myBaits", "Twist"))) 
df1 <- data.frame(Frequency = shotgunAD_filtered%>%filter(CG>0)%>%select(CG), Group = "shotgun")
df2 <- data.frame(Frequency = mybaitsAD_filtered%>%filter(CG>0)%>%select(CG), Group = "myBaits")
df3 <- data.frame(Frequency = twistAD_filtered%>%filter(CG>0)%>%select(CG), Group = "Twist")
combined_CG <- rbind(df1, df2, df3) %>% rename("Frequency"=1) %>%
  mutate(Group = factor(Group, levels = c("shotgun", "myBaits", "Twist"))) 
df1 <- data.frame(Frequency = shotgunAD_filtered%>%filter(CT>0)%>%select(CT), Group = "shotgun")
df2 <- data.frame(Frequency = mybaitsAD_filtered%>%filter(CT>0)%>%select(CT), Group = "myBaits")
df3 <- data.frame(Frequency = twistAD_filtered%>%filter(CT>0)%>%select(CT), Group = "Twist")
combined_CT <- rbind(df1, df2, df3) %>% rename("Frequency"=1) %>%
  mutate(Group = factor(Group, levels = c("shotgun", "myBaits", "Twist"))) 
df1 <- data.frame(Frequency = shotgunAD_filtered%>%filter(GT>0)%>%select(GT), Group = "shotgun")
df2 <- data.frame(Frequency = mybaitsAD_filtered%>%filter(GT>0)%>%select(GT), Group = "myBaits")
df3 <- data.frame(Frequency = twistAD_filtered%>%filter(GT>0)%>%select(GT), Group = "Twist")
combined_GT <- rbind(df1, df2, df3) %>% rename("Frequency"=1) %>%
  mutate(Group = factor(Group, levels = c("shotgun", "myBaits", "Twist"))) 
c1 <- ggplot(combined_AG, aes(x = Group, y = Frequency, fill = Group)) +
  geom_violin(trim = TRUE) +  # Use trim = TRUE if you want to trim the tails
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +  # Add boxplot inside violin
  scale_fill_manual(values = c("orange","#F8766D", "#00BFC4"), name = "Method", labels = c("shotgun","myBaits", "Twist")) +
  labs(y="NCP002") +
  theme_minimal() +
  theme(legend.position = "none") + 
  theme(axis.title.x = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),axis.line.x = element_blank())
c2 <- ggplot(combined_CT, aes(x = Group, y = Frequency, fill = Group)) +
  geom_violin(trim = TRUE) +  # Use trim = TRUE if you want to trim the tails
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +  
  scale_fill_manual(values = c("orange","#F8766D", "#00BFC4"), name = "Method", labels = c("shotgun","myBaits", "Twist")) +
  theme_minimal() +
  theme(legend.position = "none") + 
  theme(axis.title.x = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),axis.line.x = element_blank(),axis.title.y = element_blank())
c3 <- ggplot(combined_AC, aes(x = Group, y = Frequency, fill = Group)) +
  geom_violin(trim = TRUE) +  # Use trim = TRUE if you want to trim the tails
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +  
  scale_fill_manual(values = c("orange","#F8766D", "#00BFC4"), name = "Method", labels = c("shotgun","myBaits", "Twist")) +
  theme_minimal() +
  theme(legend.position = "none") + 
  theme(axis.title.x = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),axis.line.x = element_blank(),axis.title.y = element_blank())
c4 <- ggplot(combined_AT, aes(x = Group, y = Frequency, fill = Group)) +
  geom_violin(trim = TRUE) +  # Use trim = TRUE if you want to trim the tails
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +  
  scale_fill_manual(values = c("orange","#F8766D", "#00BFC4"), name = "Method", labels = c("shotgun","myBaits", "Twist")) +
  theme_minimal() +
  theme(legend.position = "none") + 
  theme(axis.title.x = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),axis.line.x = element_blank(),axis.title.y = element_blank())
c5 <- ggplot(combined_CG, aes(x = Group, y = Frequency, fill = Group)) +
  geom_violin(trim = TRUE) +  # Use trim = TRUE if you want to trim the tails
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +  
  scale_fill_manual(values = c("orange","#F8766D", "#00BFC4"), name = "Method", labels = c("shotgun","myBaits", "Twist")) +
  theme_minimal() +
  theme(legend.position = "none") + 
  theme(axis.title.x = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),axis.line.x = element_blank(),axis.title.y = element_blank())
c6 <- ggplot(combined_GT, aes(x = Group, y = Frequency, fill = Group)) +
  geom_violin(trim = TRUE) +  # Use trim = TRUE if you want to trim the tails
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +  
  scale_fill_manual(values = c("orange","#F8766D", "#00BFC4"), name = "Method", labels = c("shotgun","myBaits", "Twist")) +
  theme_minimal() +
  theme(legend.position = "none") + 
  theme(axis.title.x = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),axis.line.x = element_blank(),axis.title.y = element_blank())

#Sample VKP001
mybaitsAD <- cbind(mybaitspos, mybaitscount[13:16]) %>%
  rename(A=4, C=5, G=6, T=7) %>%
  mutate(totDepth=A+C+G+T)  %>%
  mutate(scaff = paste(chr, pos, sep = "_")) %>%
  filter(totDepth >2) %>%
  rowwise() %>%
  filter(sum(c(A > 0, C > 0, G > 0, T > 0)) == 2) %>%
  mutate(
    AC = ifelse(A > 0 & C > 0, A/(A+C), 0),
    AG = ifelse(A > 0 & G > 0, A/(A+G), 0),
    AT = ifelse(A > 0 & T > 0, A/(A+T), 0),
    CG = ifelse(C > 0 & G > 0, C/(C+G), 0),
    CT = ifelse(C > 0 & T > 0, C/(C+T), 0),
    GT = ifelse(G > 0 & T > 0, G/(G+T), 0))
twistAD <- cbind(twistpos, twistcount[13:16]) %>%
  rename(A=4, C=5, G=6, T=7) %>%
  mutate(totDepth=A+C+G+T)  %>%
  mutate(scaff = paste(chr, pos, sep = "_")) %>%
  filter(totDepth >2) %>%
  rowwise() %>%
  filter(sum(c(A > 0, C > 0, G > 0, T > 0)) == 2) %>%
  mutate(
    AC = ifelse(A > 0 & C > 0, A/(A+C), 0),
    AG = ifelse(A > 0 & G > 0, A/(A+G), 0),
    AT = ifelse(A > 0 & T > 0, A/(A+T), 0),
    CG = ifelse(C > 0 & G > 0, C/(C+G), 0),
    CT = ifelse(C > 0 & T > 0, C/(C+T), 0),
    GT = ifelse(G > 0 & T > 0, G/(G+T), 0))
shotgunAD <- cbind(shotgunpos, shotguncount[13:16]) %>%
  rename(A=4, C=5, G=6, T=7) %>%
  mutate(totDepth=A+C+G+T)  %>%
  mutate(scaff = paste(chr, pos, sep = "_")) %>%
  filter(totDepth >2) %>%
  rowwise() %>%
  filter(sum(c(A > 0, C > 0, G > 0, T > 0)) == 2) %>%
  mutate(
    AC = ifelse(A > 0 & C > 0, A/(A+C), 0),
    AG = ifelse(A > 0 & G > 0, A/(A+G), 0),
    AT = ifelse(A > 0 & T > 0, A/(A+T), 0),
    CG = ifelse(C > 0 & G > 0, C/(C+G), 0),
    CT = ifelse(C > 0 & T > 0, C/(C+T), 0),
    GT = ifelse(G > 0 & T > 0, G/(G+T), 0))
common_rows <- Reduce(intersect, list(rownames(het_vkp001), mybaitsAD$scaff, twistAD$scaff, shotgunAD$scaff)) #8879
mybaitsAD_filtered <- mybaitsAD[mybaitsAD$scaff %in% common_rows, ]
twistAD_filtered <- twistAD[twistAD$scaff %in% common_rows, ]
shotgunAD_filtered <- shotgunAD[shotgunAD$scaff %in% common_rows, ]

res_mybaits <- data.frame(Genotype=character(), estADratio=numeric(), p.value=numeric(), 
                          CI.low=numeric(), CI.high=numeric(), stringsAsFactors=FALSE)
logistic<-function(y=0){1/(1+exp(-y))}
for (i in c("AC", "AG", "AT", "CG", "CT", "GT")) {
  allele1 <- substr(i,1,1)
  allele2 <- substr(i,2,2)
  filtered_data <- mybaitsAD_filtered %>% filter(.data[[i]] > 0)
  allele1_counts <- filtered_data[[allele1]]
  allele2_counts <- filtered_data[[allele2]]
  lm1 <- glm(cbind(allele1_counts, allele2_counts)~1,data=filtered_data,family="quasibinomial")
  p.value<-summary(lm1)$coefficients[4]
  p.trans<-logistic(lm1$coefficients)
  CIlow.trans<-logistic(lm1$coefficients-summary(lm1)$coefficients[2]*1.96)
  CIup.trans<-logistic(lm1$coefficients+summary(lm1)$coefficients[2]*1.96)
  res_mybaits <- rbind(res_mybaits, data.frame(Genotype=i, estADratio=p.trans, p.value=p.value, 
                                               CI.low=CIlow.trans, CI.high=CIup.trans)) }

res_twist <- data.frame(Genotype=character(), estADratio=numeric(), p.value=numeric(), 
                        CI.low=numeric(), CI.high=numeric(), stringsAsFactors=FALSE)
logistic<-function(y=0){1/(1+exp(-y))}
for (i in c("AC", "AG", "AT", "CG", "CT", "GT")) {
  allele1 <- substr(i,1,1)
  allele2 <- substr(i,2,2)
  filtered_data <- twistAD_filtered %>% filter(.data[[i]] > 0)
  allele1_counts <- filtered_data[[allele1]]
  allele2_counts <- filtered_data[[allele2]]
  lm1 <- glm(cbind(allele1_counts, allele2_counts)~1,data=filtered_data,family="quasibinomial")
  p.value<-summary(lm1)$coefficients[4]
  p.trans<-logistic(lm1$coefficients)
  CIlow.trans<-logistic(lm1$coefficients-summary(lm1)$coefficients[2]*1.96)
  CIup.trans<-logistic(lm1$coefficients+summary(lm1)$coefficients[2]*1.96)
  res_twist <- rbind(res_twist, data.frame(Genotype=i, estADratio=p.trans, p.value=p.value, 
                                           CI.low=CIlow.trans, CI.high=CIup.trans)) }

res_shotgun <- data.frame(Genotype=character(), estADratio=numeric(), p.value=numeric(), 
                          CI.low=numeric(), CI.high=numeric(), stringsAsFactors=FALSE)
logistic<-function(y=0){1/(1+exp(-y))}
for (i in c("AC", "AG", "AT", "CG", "CT", "GT")) {
  allele1 <- substr(i,1,1)
  allele2 <- substr(i,2,2)
  filtered_data <- shotgunAD_filtered %>% filter(.data[[i]] > 0)
  allele1_counts <- filtered_data[[allele1]]
  allele2_counts <- filtered_data[[allele2]]
  lm1 <- glm(cbind(allele1_counts, allele2_counts)~1,data=filtered_data,family="quasibinomial")
  p.value<-summary(lm1)$coefficients[4]
  p.trans<-logistic(lm1$coefficients)
  CIlow.trans<-logistic(lm1$coefficients-summary(lm1)$coefficients[2]*1.96)
  CIup.trans<-logistic(lm1$coefficients+summary(lm1)$coefficients[2]*1.96)
  res_shotgun <- rbind(res_shotgun, data.frame(Genotype=i, estADratio=p.trans, p.value=p.value, 
                                               CI.low=CIlow.trans, CI.high=CIup.trans)) }

res_shotgun
res_mybaits
res_twist
#### VIOLIN PLOT OF ALLELE BIAS
df1 <- data.frame(Frequency = shotgunAD_filtered%>%filter(AC>0)%>%select(AC), Group = "shotgun")
df2 <- data.frame(Frequency = mybaitsAD_filtered%>%filter(AC>0)%>%select(AC), Group = "myBaits")
df3 <- data.frame(Frequency = twistAD_filtered%>%filter(AC>0)%>%select(AC), Group = "Twist")
combined_AC <- rbind(df1, df2, df3) %>% rename("Frequency"=1) %>%
  mutate(Group = factor(Group, levels = c("shotgun", "myBaits", "Twist"))) 
df1 <- data.frame(Frequency = shotgunAD_filtered%>%filter(AG>0)%>%select(AG), Group = "shotgun")
df2 <- data.frame(Frequency = mybaitsAD_filtered%>%filter(AG>0)%>%select(AG), Group = "myBaits")
df3 <- data.frame(Frequency = twistAD_filtered%>%filter(AG>0)%>%select(AG), Group = "Twist")
combined_AG <- rbind(df1, df2, df3) %>% rename("Frequency"=1) %>%
  mutate(Group = factor(Group, levels = c("shotgun", "myBaits", "Twist"))) 
df1 <- data.frame(Frequency = shotgunAD_filtered%>%filter(AT>0)%>%select(AT), Group = "shotgun")
df2 <- data.frame(Frequency = mybaitsAD_filtered%>%filter(AT>0)%>%select(AT), Group = "myBaits")
df3 <- data.frame(Frequency = twistAD_filtered%>%filter(AT>0)%>%select(AT), Group = "Twist")
combined_AT <- rbind(df1, df2, df3) %>% rename("Frequency"=1) %>%
  mutate(Group = factor(Group, levels = c("shotgun", "myBaits", "Twist"))) 
df1 <- data.frame(Frequency = shotgunAD_filtered%>%filter(CG>0)%>%select(CG), Group = "shotgun")
df2 <- data.frame(Frequency = mybaitsAD_filtered%>%filter(CG>0)%>%select(CG), Group = "myBaits")
df3 <- data.frame(Frequency = twistAD_filtered%>%filter(CG>0)%>%select(CG), Group = "Twist")
combined_CG <- rbind(df1, df2, df3) %>% rename("Frequency"=1) %>%
  mutate(Group = factor(Group, levels = c("shotgun", "myBaits", "Twist"))) 
df1 <- data.frame(Frequency = shotgunAD_filtered%>%filter(CT>0)%>%select(CT), Group = "shotgun")
df2 <- data.frame(Frequency = mybaitsAD_filtered%>%filter(CT>0)%>%select(CT), Group = "myBaits")
df3 <- data.frame(Frequency = twistAD_filtered%>%filter(CT>0)%>%select(CT), Group = "Twist")
combined_CT <- rbind(df1, df2, df3) %>% rename("Frequency"=1) %>%
  mutate(Group = factor(Group, levels = c("shotgun", "myBaits", "Twist"))) 
df1 <- data.frame(Frequency = shotgunAD_filtered%>%filter(GT>0)%>%select(GT), Group = "shotgun")
df2 <- data.frame(Frequency = mybaitsAD_filtered%>%filter(GT>0)%>%select(GT), Group = "myBaits")
df3 <- data.frame(Frequency = twistAD_filtered%>%filter(GT>0)%>%select(GT), Group = "Twist")
combined_GT <- rbind(df1, df2, df3) %>% rename("Frequency"=1) %>%
  mutate(Group = factor(Group, levels = c("shotgun", "myBaits", "Twist"))) 
d1 <- ggplot(combined_AG, aes(x = Group, y = Frequency, fill = Group)) +
  geom_violin(trim = TRUE) +  # Use trim = TRUE if you want to trim the tails
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +  # Add boxplot inside violin
  scale_fill_manual(values = c("orange","#F8766D", "#00BFC4"), name = "Method", labels = c("Shotgun","myBaits", "Twist")) +
  labs(y="VKP001") +
  theme_minimal() +
  theme(legend.position = "none") + 
  theme(axis.title.x = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),axis.line.x = element_blank())
d2 <- ggplot(combined_CT, aes(x = Group, y = Frequency, fill = Group)) +
  geom_violin(trim = TRUE) +  # Use trim = TRUE if you want to trim the tails
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +  
  scale_fill_manual(values = c("orange","#F8766D", "#00BFC4"), name = "Method", labels = c("Shotgun","myBaits", "Twist")) +
  theme_minimal() +
  theme(legend.position = "none") + 
  theme(axis.title.x = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),axis.line.x = element_blank(),axis.title.y = element_blank())
d3 <- ggplot(combined_AC, aes(x = Group, y = Frequency, fill = Group)) +
  geom_violin(trim = TRUE) +  # Use trim = TRUE if you want to trim the tails
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +  
  scale_fill_manual(values = c("orange","#F8766D", "#00BFC4"), name = "Method", labels = c("Shotgun","myBaits", "Twist")) +
  theme_minimal() +
  theme(legend.position = "none") + 
  theme(axis.title.x = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),axis.line.x = element_blank(),axis.title.y = element_blank())
d4 <- ggplot(combined_AT, aes(x = Group, y = Frequency, fill = Group)) +
  geom_violin(trim = TRUE) +  # Use trim = TRUE if you want to trim the tails
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +  
  scale_fill_manual(values = c("orange","#F8766D", "#00BFC4"), name = "Method", labels = c("Shotgun","myBaits", "Twist")) +
  theme_minimal() +
  theme(legend.position = "none") + 
  theme(axis.title.x = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),axis.line.x = element_blank(),axis.title.y = element_blank())
d5 <- ggplot(combined_CG, aes(x = Group, y = Frequency, fill = Group)) +
  geom_violin(trim = TRUE) +  # Use trim = TRUE if you want to trim the tails
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +  
  scale_fill_manual(values = c("orange","#F8766D", "#00BFC4"), name = "Method", labels = c("Shotgun","myBaits", "Twist")) +
  theme_minimal() +
  theme(legend.position = "none") + 
  theme(axis.title.x = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),axis.line.x = element_blank(),axis.title.y = element_blank())
d6 <- ggplot(combined_GT, aes(x = Group, y = Frequency, fill = Group)) +
  geom_violin(trim = TRUE) +  # Use trim = TRUE if you want to trim the tails
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +  
  scale_fill_manual(values = c("orange","#F8766D", "#00BFC4"), name = "Method", labels = c("Shotgun","myBaits", "Twist")) +
  theme_minimal() +
  theme(legend.position = "none") + 
  theme(axis.title.x = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),axis.line.x = element_blank(),axis.title.y = element_blank())

# Create titles for each column
title1 <- textGrob("AG", gp = gpar(fontsize = 14, fontface = "bold"))
title2 <- textGrob("CT", gp = gpar(fontsize = 14, fontface = "bold"))
title3 <- textGrob("AC", gp = gpar(fontsize = 14, fontface = "bold"))
title4 <- textGrob("AT", gp = gpar(fontsize = 14, fontface = "bold"))
title5 <- textGrob("CG", gp = gpar(fontsize = 14, fontface = "bold"))
title6 <- textGrob("GT", gp = gpar(fontsize = 14, fontface = "bold"))

grob_plots <- arrangeGrob(a1,a2,a3,a4,a5,a6,
                          d1,d2,d3,d4,d5,d6,
                          b1,b2,b3,b4,b5,b6,
                          c1,c2,c3,c4,c5,c6, nrow=4,ncol=6,widths = c(1,1,1,1,1,1))

combined <- grid.arrange(arrangeGrob(title1, title2, title3, title4, title5, title6, ncol = 6, widths = c(1,1,1,1,1,1)),
                         grob_plots, shared_legend, nrow = 3, heights = c(1,10,1))
#14x7
