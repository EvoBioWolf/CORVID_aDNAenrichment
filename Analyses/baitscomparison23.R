library(ggplot2)
library(readxl)
library(dplyr)
library(gridExtra)
library(grid)
library(tibble)
library(ggvenn)
library(data.table)
library(tidyr)
library(ggsignif)
library(patchwork)
library(ggrepel)

#setwd("/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/05_aDNA/00_baitscomparison") 
setwd("~/Dropbox/CORVID_baits/Analyses") 
dat <- read.csv("baitscomparison23.csv",header=TRUE, sep=",") %>%
  mutate(All_mapped_reads_to_original_panel = ifelse(Method == "myBaits", All_mapped_reads_to_104k_121bp, All_mapped_reads_to_232k_80bp)) %>%
  mutate(Nr_comparable_mapped_reads = Nr_mapped_reads - All_non_comparable_mapped_reads) %>%
  mutate(Target_eff_adjusted = All_mapped_reads_to_104k_80bp / Nr_comparable_mapped_reads*100) %>%
  mutate(Target_eff_original = All_mapped_reads_to_original_panel / Nr_mapped_reads*100) %>%
  mutate(Target_eff_adjusted_more = ifelse(Method == "Twist", (All_mapped_reads_to_104k_80bp / (Nr_comparable_mapped_reads*104/232))*100, Target_eff_adjusted)) %>%
  mutate(Target_eff_adjusted_more2 = ifelse(Method == "Twist", (All_mapped_reads_to_104k_80bp*2 / (Nr_comparable_mapped_reads*104/232))*100, Target_eff_adjusted)) %>%
  mutate(Target_eff_original_2x = ifelse(Method == "Twist", (All_mapped_reads_to_original_panel*2 / Nr_mapped_reads*100), Target_eff_original)) %>%
  mutate(Ontarget_rate_adjusted = All_mapped_reads_to_104k_80bp / (Nr_rawreads - All_non_comparable_mapped_reads)*100) %>%
  mutate(Ontarget_rate = All_mapped_reads_to_original_panel / Nr_rawreads*100) %>%
  mutate(Obs_exp_cov = Mean_cov_of_104k_SNPsite_incdup / Expected_genomic_coverage_of_input) %>%
  mutate(Percentage_of_mtDNA = All_MT_reads / Nr_mapped_reads * 100)

# Figure 2: Capture efficiency and fold enrichment --------------------------------------------------------------
scatter_data <- dat %>%
  select(Sample, Method, Endogenous_DNA, Target_eff_adjusted, Target_eff_original, Ontarget_rate_adjusted, Ontarget_rate, 
         Target_eff_adjusted_more, Target_eff_adjusted_more2,Target_eff_original_2x,Percentage_of_mtDNA,
         Obs_exp_cov, Mean_cov_of_104k_SNPsite_incdup, Mean_cov_of_target_104k_SNPsite_dedup,Mean_cov_of_104k_80bp_incdup) %>%
  pivot_wider(names_from = Method, values_from = c(Endogenous_DNA,Target_eff_adjusted, Target_eff_original, Ontarget_rate_adjusted, Ontarget_rate, 
                                                   Target_eff_adjusted_more, Target_eff_adjusted_more2,Target_eff_original_2x,Percentage_of_mtDNA,
                                                   Obs_exp_cov, Mean_cov_of_104k_SNPsite_incdup, Mean_cov_of_target_104k_SNPsite_dedup,Mean_cov_of_104k_80bp_incdup))

model <- lm(Endogenous_DNA_Twist ~ Endogenous_DNA_myBaits, data = scatter_data)
# Test if slope is different from 1
summary_model <- summary(model)
slope <- coef(model)[2]
se_slope <- summary_model$coefficients[2, 2]
t_slope <- (slope - 1) / se_slope
p_value_slope_vs_1 <- 2 * pt(-abs(t_slope), df = model$df.residual)
# Test if intercept is different from 0
intercept <- coef(model)[1]
se_intercept <- summary_model$coefficients[1, 2]
t_intercept <- (intercept - 0) / se_intercept
p_value_intercept_vs_0 <- 2 * pt(-abs(t_intercept), df = model$df.residual)
p_value <- summary(model)$coefficients[2, 4]
r2 <- summary(model)$r.squared
eqn <- paste0("\nR² = ", round(r2, 3),
  "\np (slope0) = ", format.pval(p_value, digits = 3, eps = 0.001),
  "\np (slope1) = ", format.pval(p_value_slope_vs_1, digits = 3, eps = 0.001),
  "\np (int0) = ", format.pval(p_value_intercept_vs_0, digits = 3, eps = 0.001))
p1 <- ggplot(scatter_data, aes(x=Endogenous_DNA_Twist , y =Endogenous_DNA_myBaits)) +
  geom_point(size = 1) + geom_smooth(method = "lm", se = TRUE, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  annotate("text", x = 5, y = 90, label = eqn, hjust = 0, size = 2.5) +
  labs(x = "Twist", y = "myBaits", title = "Reads mapped to\n crow genome (%)") + ylim(0,100) + xlim(0,100) +
  theme_minimal() +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8))

model <- lm(Target_eff_adjusted_Twist ~ Target_eff_adjusted_myBaits, data = scatter_data)
summary_model <- summary(model)
slope <- coef(model)[2]
se_slope <- summary_model$coefficients[2, 2]
t_slope <- (slope - 1) / se_slope
p_value_slope_vs_1 <- 2 * pt(-abs(t_slope), df = model$df.residual)
intercept <- coef(model)[1]
se_intercept <- summary_model$coefficients[1, 2]
t_intercept <- (intercept - 0) / se_intercept
p_value_intercept_vs_0 <- 2 * pt(-abs(t_intercept), df = model$df.residual)
p_value <- summary(model)$coefficients[2, 4]
r2 <- summary(model)$r.squared
eqn <- paste0("\nR² = ", round(r2, 3),
              "\np (slope0) = ", format.pval(p_value, digits = 3, eps = 0.001),
              "\np (slope1) = ", format.pval(p_value_slope_vs_1, digits = 3, eps = 0.001),
              "\np (int0) = ", format.pval(p_value_intercept_vs_0, digits = 3, eps = 0.001))
p2 <- ggplot(scatter_data, aes(x=Target_eff_adjusted_Twist , y =Target_eff_adjusted_myBaits)) +
  geom_point(size = 1) + geom_smooth(method = "lm", se = TRUE, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  annotate("text", x = 30, y = 25, label = eqn, hjust = 0, size = 2.5) +
  labs(x = "Twist", y = "myBaits", title = "Target efficiency (%)\n (comparable panel)") + ylim(0,100) + xlim(0,100) +
  theme_minimal() +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8))

model <- lm(Target_eff_adjusted_more_Twist ~ Target_eff_adjusted_more_myBaits, data = scatter_data)
summary_model <- summary(model)
slope <- coef(model)[2]
se_slope <- summary_model$coefficients[2, 2]
t_slope <- (slope - 1) / se_slope
p_value_slope_vs_1 <- 2 * pt(-abs(t_slope), df = model$df.residual)
intercept <- coef(model)[1]
se_intercept <- summary_model$coefficients[1, 2]
t_intercept <- (intercept - 0) / se_intercept
p_value_intercept_vs_0 <- 2 * pt(-abs(t_intercept), df = model$df.residual)
p_value <- summary(model)$coefficients[2, 4]
r2 <- summary(model)$r.squared
eqn <- paste0("\nR² = ", round(r2, 3),
              "\np (slope0) = ", format.pval(p_value, digits = 3, eps = 0.001),
              "\np (slope1) = ", format.pval(p_value_slope_vs_1, digits = 3, eps = 0.001),
              "\np (int0) = ", format.pval(p_value_intercept_vs_0, digits = 3, eps = 0.001))
supp_x <- ggplot(scatter_data, aes(x=Target_eff_adjusted_more_Twist , y =Target_eff_adjusted_more_myBaits)) +
  geom_point(size = 1) + geom_smooth(method = "lm", se = TRUE, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  annotate("text", x = 50, y = 30, label = eqn, hjust = 0, size = 2.5) +
  labs(x = "Twist", y = "myBaits", title = "Target efficiency of comparable panel (%)\n (104/232 off-target for Twist)") + ylim(0,100) + xlim(0,100) +
  theme_minimal() +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8))

model <- lm(Target_eff_adjusted_more2_Twist ~ Target_eff_adjusted_more2_myBaits, data = scatter_data)
summary_model <- summary(model)
slope <- coef(model)[2]
se_slope <- summary_model$coefficients[2, 2]
t_slope <- (slope - 1) / se_slope
p_value_slope_vs_1 <- 2 * pt(-abs(t_slope), df = model$df.residual)
intercept <- coef(model)[1]
se_intercept <- summary_model$coefficients[1, 2]
t_intercept <- (intercept - 0) / se_intercept
p_value_intercept_vs_0 <- 2 * pt(-abs(t_intercept), df = model$df.residual)
p_value <- summary(model)$coefficients[2, 4]
r2 <- summary(model)$r.squared
eqn <- paste0("\nR² = ", round(r2, 3),
              "\np (slope0) = ", format.pval(p_value, digits = 3, eps = 0.001),
              "\np (slope1) = ", format.pval(p_value_slope_vs_1, digits = 3, eps = 0.001),
              "\np (int0) = ", format.pval(p_value_intercept_vs_0, digits = 3, eps = 0.001))
supp_y <- ggplot(scatter_data, aes(x=Target_eff_adjusted_more2_Twist , y =Target_eff_adjusted_more2_myBaits)) +
  geom_point(size = 1) + geom_smooth(method = "lm", se = TRUE, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  annotate("text", x = 50, y = 30, label = eqn, hjust = 0, size = 2.5) +
  labs(x = "Twist", y = "myBaits", title = "Target efficiency of comparable panel (%)\n (2X on-target & 104/232 off-target for Twist)") + ylim(0,100) + xlim(0,100) +
  theme_minimal() +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8))

model <- lm(Target_eff_original_2x_Twist ~ Target_eff_original_2x_myBaits, data = scatter_data)
summary_model <- summary(model)
slope <- coef(model)[2]
se_slope <- summary_model$coefficients[2, 2]
t_slope <- (slope - 1) / se_slope
p_value_slope_vs_1 <- 2 * pt(-abs(t_slope), df = model$df.residual)
intercept <- coef(model)[1]
se_intercept <- summary_model$coefficients[1, 2]
t_intercept <- (intercept - 0) / se_intercept
p_value_intercept_vs_0 <- 2 * pt(-abs(t_intercept), df = model$df.residual)
p_value <- summary(model)$coefficients[2, 4]
r2 <- summary(model)$r.squared
eqn <- paste0("\nR² = ", round(r2, 3),
              "\np (slope0) = ", format.pval(p_value, digits = 3, eps = 0.001),
              "\np (slope1) = ", format.pval(p_value_slope_vs_1, digits = 3, eps = 0.001),
              "\np (int0) = ", format.pval(p_value_intercept_vs_0, digits = 3, eps = 0.001))
supp_z <- ggplot(scatter_data, aes(x=Target_eff_original_2x_Twist , y =Target_eff_original_2x_myBaits)) +
  geom_point(size = 1) + geom_smooth(method = "lm", se = TRUE, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  annotate("text", x = 50, y = 30, label = eqn, hjust = 0, size = 2.5) +
  labs(x = "Twist", y = "myBaits", title = "Target efficiency of original panel (%)\n (2X on-target for Twist)") + ylim(0,100) + xlim(0,100) +
  theme_minimal() +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8))

model <- lm(Target_eff_original_Twist ~ Target_eff_original_myBaits, data = scatter_data)
summary_model <- summary(model)
slope <- coef(model)[2]
se_slope <- summary_model$coefficients[2, 2]
t_slope <- (slope - 1) / se_slope
p_value_slope_vs_1 <- 2 * pt(-abs(t_slope), df = model$df.residual)
intercept <- coef(model)[1]
se_intercept <- summary_model$coefficients[1, 2]
t_intercept <- (intercept - 0) / se_intercept
p_value_intercept_vs_0 <- 2 * pt(-abs(t_intercept), df = model$df.residual)
p_value <- summary(model)$coefficients[2, 4]
r2 <- summary(model)$r.squared
eqn <- paste0("\nR² = ", round(r2, 3),
              "\np (slope0) = ", format.pval(p_value, digits = 3, eps = 0.001),
              "\np (slope1) = ", format.pval(p_value_slope_vs_1, digits = 3, eps = 0.001),
              "\np (int0) = ", format.pval(p_value_intercept_vs_0, digits = 3, eps = 0.001))
p3 <- ggplot(scatter_data, aes(x=Target_eff_original_Twist , y =Target_eff_original_myBaits)) +
  geom_point(size = 1) + geom_smooth(method = "lm", se = TRUE, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  annotate("text", x = 30, y = 25, label = eqn, hjust = 0, size = 2.5) +
  labs(x = "Twist", y = "myBaits", title = "Target efficiency (%)\n (original panel)") + ylim(0,100) + xlim(0,100) +
  theme_minimal() +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8))

model <- lm(Ontarget_rate_Twist ~ Ontarget_rate_myBaits, data = scatter_data)
summary_model <- summary(model)
slope <- coef(model)[2]
se_slope <- summary_model$coefficients[2, 2]
t_slope <- (slope - 1) / se_slope
p_value_slope_vs_1 <- 2 * pt(-abs(t_slope), df = model$df.residual)
intercept <- coef(model)[1]
se_intercept <- summary_model$coefficients[1, 2]
t_intercept <- (intercept - 0) / se_intercept
p_value_intercept_vs_0 <- 2 * pt(-abs(t_intercept), df = model$df.residual)
p_value <- summary(model)$coefficients[2, 4]
r2 <- summary(model)$r.squared
eqn <- paste0("\nR² = ", round(r2, 3),
              "\np (slope0) = ", format.pval(p_value, digits = 3, eps = 0.001),
              "\np (slope1) = ", format.pval(p_value_slope_vs_1, digits = 3, eps = 0.001),
              "\np (int0) = ", format.pval(p_value_intercept_vs_0, digits = 3, eps = 0.001))
p4 <- ggplot(scatter_data, aes(x=Ontarget_rate_Twist , y =Ontarget_rate_myBaits)) +
  geom_point(size = 1) + geom_smooth(method = "lm", se = TRUE, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  annotate("text", x = 5, y = 90, label = eqn, hjust = 0, size = 2.5) +
  labs(x = "Twist", y = "myBaits", title = "On-target rate (%)\n (original panel)") + ylim(0,100) + xlim(0,100) + 
  theme_minimal() +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8))

model <- lm(Percentage_of_mtDNA_Twist ~ Percentage_of_mtDNA_myBaits, data = scatter_data)
summary_model <- summary(model)
slope <- coef(model)[2]
se_slope <- summary_model$coefficients[2, 2]
t_slope <- (slope - 1) / se_slope
p_value_slope_vs_1 <- 2 * pt(-abs(t_slope), df = model$df.residual)
intercept <- coef(model)[1]
se_intercept <- summary_model$coefficients[1, 2]
t_intercept <- (intercept - 0) / se_intercept
p_value_intercept_vs_0 <- 2 * pt(-abs(t_intercept), df = model$df.residual)
p_value <- summary(model)$coefficients[2, 4]
r2 <- summary(model)$r.squared
eqn <- paste0("\nR² = ", round(r2, 3),
              "\np (slope0) = ", format.pval(p_value, digits = 3, eps = 0.001),
              "\np (slope1) = ", format.pval(p_value_slope_vs_1, digits = 3, eps = 0.001),
              "\np (int0) = ", format.pval(p_value_intercept_vs_0, digits = 3, eps = 0.001))
p5 <- ggplot(scatter_data, aes(x=Percentage_of_mtDNA_Twist , y =Percentage_of_mtDNA_myBaits)) +
  geom_point(size = 1) + geom_smooth(method = "lm", se = TRUE, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  annotate("text", x = 0.05, y = 0.03, label = eqn, hjust = 0, size = 2.5) +
  geom_text(data = scatter_data %>% filter(Percentage_of_mtDNA_myBaits>0.05),aes(label = Sample),size = 2.5, hjust=0, vjust = -1) +
  labs(x = "Twist", y = "myBaits", title = "Mito reads among all\n mapped reads (%)") + ylim(0,0.11) + xlim(0,0.11) +
  theme_minimal() +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8))

Wgs_mapped_to_104k_80bp <- rep(c(20820173,1927648,307075,2811357),2)
Wgs_mapped_to_104k_121bp <- rep(c(27475838,2590139,408977,3716873),2)
Wgs_mapped_to_232k_80bp <- rep(c(47585894,4437366,691157,6387808),2)
WGS_all_mapped_reads <- rep(c(1642132057, 161971725, 24544200, 223432035),2)
Wgs_non_comparable_wgs_reads = (Wgs_mapped_to_104k_121bp-Wgs_mapped_to_104k_80bp)+(Wgs_mapped_to_232k_80bp-Wgs_mapped_to_104k_80bp)

fold_eff_dat <- dat %>% filter(Sample %in% c("BRW001","DVT014","NCP002","VKP001")) %>% 
  mutate(Wgs_mapped_to_104k_80bp = Wgs_mapped_to_104k_80bp,
         Wgs_mapped_to_104k_121bp = Wgs_mapped_to_104k_121bp,
         Wgs_mapped_to_232k_80bp = Wgs_mapped_to_232k_80bp,
         Wgs_non_comparable_wgs_reads = Wgs_non_comparable_wgs_reads) %>%
  mutate(Wgs_target_eff_adjusted = Wgs_mapped_to_104k_80bp/(WGS_all_mapped_reads-Wgs_non_comparable_wgs_reads)*100) %>%
  mutate(Wgs_target_eff_original =  ifelse(Method == "myBaits", Wgs_mapped_to_104k_121bp/WGS_all_mapped_reads*100, Wgs_mapped_to_232k_80bp/WGS_all_mapped_reads*100)) %>%
  mutate(fold_enrich_adjusted = Target_eff_adjusted / Wgs_target_eff_adjusted) %>%
  mutate(fold_enrich_original = Target_eff_original / Wgs_target_eff_original)

fold1 <- ggplot(fold_eff_dat, aes(x = Sample, y = fold_enrich_adjusted, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Sample", y = "Fold Enrichment", title="Comparable panel") +
  theme_minimal() + ylim(0, 60) +  guides(fill = guide_legend(ncol = 2)) +
  scale_fill_manual(values = c("myBaits" = "#E69F00", "Twist" = "#0072B2")) + theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),legend.position = c(0.5, 1),legend.title = element_blank())

fold2 <- ggplot(fold_eff_dat, aes(x = Sample, y = fold_enrich_original, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Sample", y = "Fold Enrichment", title = "Original panel") +
  theme_minimal() + ylim(0, 60) + guides(fill = guide_legend(ncol = 2)) +
  scale_fill_manual(values = c("myBaits" = "#E69F00", "Twist" = "#0072B2")) + theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),legend.position = c(0.5, 1),legend.title = element_blank())

row1 <- p4 + p1 + p2 + p3 + plot_layout(widths = c(1, 1, 1, 1))
row2 <- fold1 + fold2 + p5 + plot_layout(widths = c(1.5, 1.5,1))
row1 / row2 + plot_annotation(tag_levels = 'A') & theme(plot.margin = unit(c(1, 1, 1, 1), "pt")) #8x5
supp <- (supp_x + supp_y) / (supp_z + plot_spacer()) + plot_annotation(tag_levels = 'A') + plot_layout(guides = 'collect', widths = c(1, 1)) #7x5

# Figure 3: Coverage ----------------------------------------------------------------------------
# Assess SNP coverage and evenness
model <- lm(Mean_cov_of_104k_SNPsite_incdup_Twist ~ Mean_cov_of_104k_SNPsite_incdup_myBaits, data = scatter_data)
summary_model <- summary(model)
slope <- coef(model)[2]
se_slope <- summary_model$coefficients[2, 2]
t_slope <- (slope - 1) / se_slope
p_value_slope_vs_1 <- 2 * pt(-abs(t_slope), df = model$df.residual)
intercept <- coef(model)[1]
se_intercept <- summary_model$coefficients[1, 2]
t_intercept <- (intercept - 0) / se_intercept
p_value_intercept_vs_0 <- 2 * pt(-abs(t_intercept), df = model$df.residual)
p_value <- summary(model)$coefficients[2, 4]
r2 <- summary(model)$r.squared
eqn <- paste0("\nR² = ", round(r2, 3),
              "\np (slope0) = ", format.pval(p_value, digits = 3, eps = 0.001),
              "\np (slope1) = ", format.pval(p_value_slope_vs_1, digits = 3, eps = 0.001),
              "\np (int0) = ", format.pval(p_value_intercept_vs_0, digits = 3, eps = 0.001))
p5 <- ggplot(scatter_data, aes(x=Mean_cov_of_104k_SNPsite_incdup_Twist , y =Mean_cov_of_104k_SNPsite_incdup_myBaits)) +
  geom_point(size = 1) + geom_smooth(method = "lm", se = TRUE, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  annotate("text", x = 100, y = 45, label = eqn, hjust = 0, size = 2.5) +
  geom_text(data = scatter_data %>% filter(Sample %in% c("BRW001", "VKP001", "TPC001")),aes(label = Sample), vjust = -1, size = 2.5) +
  labs(x = "Twist", y = "myBaits", title = "Mean coverage of 104K\n SNP sites with duplicates") + ylim(0,220) + xlim(0,220) +
  theme_minimal() +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8))

model <- lm(Mean_cov_of_target_104k_SNPsite_dedup_Twist ~ Mean_cov_of_target_104k_SNPsite_dedup_myBaits, data = scatter_data)
summary_model <- summary(model)
slope <- coef(model)[2]
se_slope <- summary_model$coefficients[2, 2]
t_slope <- (slope - 1) / se_slope
p_value_slope_vs_1 <- 2 * pt(-abs(t_slope), df = model$df.residual)
intercept <- coef(model)[1]
se_intercept <- summary_model$coefficients[1, 2]
t_intercept <- (intercept - 0) / se_intercept
p_value_intercept_vs_0 <- 2 * pt(-abs(t_intercept), df = model$df.residual)
p_value <- summary(model)$coefficients[2, 4]
r2 <- summary(model)$r.squared
eqn <- paste0("\nR² = ", round(r2, 3),
              "\np (slope0) = ", format.pval(p_value, digits = 3, eps = 0.001),
              "\np (slope1) = ", format.pval(p_value_slope_vs_1, digits = 3, eps = 0.001),
              "\np (int0) = ", format.pval(p_value_intercept_vs_0, digits = 3, eps = 0.001))
p6 <- ggplot(scatter_data, aes(x=Mean_cov_of_target_104k_SNPsite_dedup_Twist , y=Mean_cov_of_target_104k_SNPsite_dedup_myBaits)) +
  geom_point(size = 1) + geom_smooth(method = "lm", se = TRUE, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  annotate("text", x = 25, y = 10, label = eqn, hjust = 0, size = 2.5) +
  geom_text(data = scatter_data %>% filter(Sample %in% c("BRW001", "VKP001", "TPC001")),aes(label = Sample), vjust = -1, size = 2.5) +
  labs(x = "Twist", y = "myBaits", title = "Mean coverage of 104K\n SNP sites without duplicates") + ylim(0,50) + xlim(0,50) +
  theme_minimal() +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8))

model <- lm(Mean_cov_of_104k_80bp_incdup_Twist ~ Mean_cov_of_104k_80bp_incdup_myBaits, data = scatter_data)
summary_model <- summary(model)
slope <- coef(model)[2]
se_slope <- summary_model$coefficients[2, 2]
t_slope <- (slope - 1) / se_slope
p_value_slope_vs_1 <- 2 * pt(-abs(t_slope), df = model$df.residual)
intercept <- coef(model)[1]
se_intercept <- summary_model$coefficients[1, 2]
t_intercept <- (intercept - 0) / se_intercept
p_value_intercept_vs_0 <- 2 * pt(-abs(t_intercept), df = model$df.residual)
p_value <- summary(model)$coefficients[2, 4]
r2 <- summary(model)$r.squared
eqn <- paste0("\nR² = ", round(r2, 3),
              "\np (slope0) = ", format.pval(p_value, digits = 3, eps = 0.001),
              "\np (slope1) = ", format.pval(p_value_slope_vs_1, digits = 3, eps = 0.001),
              "\np (int0) = ", format.pval(p_value_intercept_vs_0, digits = 3, eps = 0.001))
p7 <- ggplot(scatter_data, aes(x=Mean_cov_of_104k_80bp_incdup_Twist , y=Mean_cov_of_104k_80bp_incdup_myBaits)) +
  geom_point(size = 1) + geom_smooth(method = "lm", se = TRUE, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  annotate("text", x = 100, y = 45, label = eqn, hjust = 0, size = 2.5) +
  geom_text(data = scatter_data %>% filter(Sample %in% c("BRW001", "VKP001", "TPC001")),aes(label = Sample), vjust = -1, size = 2.5) +
  labs(x = "Twist", y = "myBaits", title = "Mean coverage of comparable\n panel with duplicates") + ylim(0,200) + xlim(0,200) +
  theme_minimal() +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8)) #one point removed in plot (outside range) 

dat1 = dat %>% arrange(desc(Age), .by_group = TRUE)
twistdat = dat1 %>% filter(Method == "Twist") 
mybaitsdat = dat1 %>% filter(Method == "myBaits")

count_80foldcov1 <- dat1 %>%
  filter(at_least_1X_of_target_104k_SNPsite_incdup*100 > 79) %>% 
  group_by(Method) %>% 
  summarize(count_above_80 = n())
count_80foldcov2 <- dat1 %>%
  filter(at_least_2X_of_target_104k_SNPsite_incdup*100 > 79) %>% 
  group_by(Method) %>% 
  summarize(count_above_80 = n())
count_80foldcov3 <- dat1 %>%
  filter(at_least_3X_of_target_104k_SNPsite_incdup*100 > 79) %>% 
  group_by(Method) %>% 
  summarize(count_above_80 = n())
count_80foldcov4 <- dat1 %>%
  filter(at_least_4X_of_target_104k_SNPsite_incdup*100 > 79) %>% 
  group_by(Method) %>% 
  summarize(count_above_80 = n())

q1<- ggplot(dat1, aes(x = Method, y =  at_least_1X_of_target_104k_SNPsite_incdup*100, fill = Method)) +
  geom_violin(trim = TRUE, alpha = 0.6, show.legend = FALSE) +
  geom_boxplot(width=0.05, alpha = 0.8, show.legend = FALSE) +
  geom_signif(comparisons = list(c("myBaits", "Twist")),map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05),tip_length = 0,y_position=0.7) +
  geom_point(position = position_jitter(width = 0.1), size = 0.8, show.legend = FALSE) +
  labs(title="At least 1X coverage",y = "Percentage of target SNPs (%)", x="") +
  coord_cartesian(ylim = c(0, 100)) +
  geom_text(data = count_80foldcov1, aes(x=Method,y = 80, label = c(paste("N = ",count_80foldcov1[1,2],sep=""), paste("N = ",count_80foldcov1[2,2],sep=""))), 
            vjust = -0.5, hjust=1.5, size = 3, color = "black") +
  geom_text_repel(data = dat1%>%filter(Sample %in%  c("KCZ016", "KCZ004")), aes(x=Method,y = at_least_1X_of_target_104k_SNPsite_incdup*100, label = Sample),hjust=-0.1, size = 2) +
  geom_hline(yintercept = 80, linetype = "dashed", color = "black") +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 10)) + 
  scale_fill_manual(values = c("myBaits" = "#E69F00", "Twist" = "#0072B2")) +
  theme_minimal() +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8))

q2<- ggplot(dat1, aes(x = Method, y =  at_least_2X_of_target_104k_SNPsite_incdup*100, fill = Method)) +
  geom_violin(trim = TRUE, alpha = 0.6, show.legend = FALSE) +
  geom_boxplot(width=0.05, alpha = 0.8, show.legend = FALSE) +
  geom_signif(comparisons = list(c("myBaits", "Twist")),map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05),tip_length = 0,y_position=0.7) +
  geom_point(position = position_jitter(width = 0.1),size = 0.8, show.legend = FALSE) + 
  labs(title="At least 2X coverage",y = "", x="") +
  geom_text(data = count_80foldcov2, aes(x=Method,y = 80, label = c(paste("N = ",count_80foldcov2[1,2],sep=""), paste("N = ",count_80foldcov2[2,2],sep=""))), 
            vjust = -0.5, hjust=1.5, size = 3, color = "black") +
  geom_text_repel(data = dat1%>%filter(Sample %in%  c("KCZ016", "KCZ004")), aes(x=Method,y = at_least_2X_of_target_104k_SNPsite_incdup*100, label = Sample),hjust=-0.1, size = 2) +
  geom_hline(yintercept = 80, linetype = "dashed", color = "black") +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 10)) + 
  scale_fill_manual(values = c("myBaits" = "#E69F00", "Twist" = "#0072B2")) +
  theme_minimal() +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8))

q3<- ggplot(dat1, aes(x = Method, y =  at_least_3X_of_target_104k_SNPsite_incdup*100, fill = Method)) +
  geom_violin(trim = TRUE, alpha = 0.6, show.legend = FALSE) +
  geom_boxplot(width=0.05, alpha = 0.8, show.legend = FALSE) +
  geom_signif(comparisons = list(c("myBaits", "Twist")),map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05),tip_length = 0,y_position=0.7) +
  geom_point(position = position_jitter(width = 0.1),size = 0.8, show.legend = FALSE) +  
  labs(title="At least 3X coverage",y = "", x="") +
  geom_text(data = count_80foldcov3, aes(x=Method,y = 80, label = c(paste("N = ",count_80foldcov3[1,2],sep=""), paste("N = ",count_80foldcov3[2,2],sep=""))), 
            vjust = -0.5, hjust=1.5, size = 3, color = "black") +
  geom_text_repel(data = dat1%>%filter(Sample %in%  c("KCZ016", "KCZ004")), aes(x=Method,y = at_least_3X_of_target_104k_SNPsite_incdup*100, label = Sample),hjust=-0.1, size = 2) +
  geom_hline(yintercept = 80, linetype = "dashed", color = "black") +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 10)) + 
  scale_fill_manual(values = c("myBaits" = "#E69F00", "Twist" = "#0072B2")) +
  theme_minimal() +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8))

q4 <- ggplot(dat1, aes(x = Method, y =  at_least_4X_of_target_104k_SNPsite_incdup*100, fill = Method)) +
  geom_violin(trim = TRUE, alpha = 0.6, show.legend = FALSE) +
  geom_boxplot(width=0.05, alpha = 0.8, show.legend = FALSE) +
  geom_signif(comparisons = list(c("myBaits", "Twist")),map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05),tip_length = 0,y_position=0.7) +
  geom_point(position = position_jitter(width = 0.1),size = 0.8, show.legend = FALSE) +  # Add jittered points
  labs(title="At least 4X coverage",y = "", x="") +
  geom_text(data = count_80foldcov4, aes(x=Method,y = 80, label = c(paste("N = ",count_80foldcov4[1,2],sep=""), paste("N = ",count_80foldcov4[2,2],sep=""))), 
            vjust = -0.5, hjust=1.5, size = 3, color = "black") +
  geom_text_repel(data = dat1%>%filter(Sample %in%  c("KCZ016", "KCZ004")), aes(x=Method,y = at_least_4X_of_target_104k_SNPsite_incdup*100, label = Sample),hjust=-0.1, size = 2) +
  geom_hline(yintercept = 80, linetype = "dashed", color = "black") +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 10)) + 
  scale_fill_manual(values = c("myBaits" = "#E69F00", "Twist" = "#0072B2")) +
  theme_minimal() +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8))

row1 <- p7 + p5 + p6 + plot_layout(guides = 'collect', widths = c(1, 1, 1))
row2 <- q1 + q2 + q3 + plot_layout(guides = 'collect', widths = c(1,1,1))
row1 / row2 + plot_annotation(tag_levels = 'A') & theme(plot.margin = unit(c(1, 1, 1, 1), "pt"),legend.position = "bottom") #8x5

# Supplementary figure: SNP coverage distribution from samtools depth -------------------
library(purrr)
mybaitscov <- read.delim("coverage_104k_mybaits.txt",header=FALSE, sep="\t") 
twistcov <- read.delim("coverage_104k_twist.txt",header=FALSE, sep="\t") 
popbaits <- read.table("~/Dropbox/CORVID_baits/Analyses/angsd/popbaits.txt", sep="\t", header=TRUE) 
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
      fill_color <- ifelse(type == "mybaitscov", "#E69F00", "#0072B2")
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
  pdf(file=paste("coverage_distribution_104k_",type,"_selected_updated",sep=""))
  grid.arrange(grobs=plot_list[popbaits_selected$Samples], ncol=4)
  dev.off()
}

# Figure 4: Enrichment across varying GC bins ----------------------------------------------------------------------------
# Figure 4a normalized coverage across GC content ------------------------------------------------------------------------
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
  scale_color_manual(values = c("#E69F00", "#0072B2"), name = "Method", labels = c("myBaits", "Twist")) +
  labs(title = "",x = "GC content",y = "Normalized coverage") +
  theme_minimal() +
  theme(legend.position = c(0.9, 0.8),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10))

# Figure 4b-c Boxplot for GC / AT dropout ----------------------------------------------------------------------------
dropout <- read.delim("./gc/GC_AT_dropout_summary.txt",header=TRUE, sep="\t") 
dropout_selected <- dropout %>% separate(col=Sample,into=c("Sample","Method"),sep="_") %>%  
  mutate(Method=recode(Method, "TE"="Twist", "mybaits"="myBaits")) %>%
  filter(Sample %in% popbaits_selected$Samples) 

drop1<- ggplot(dropout_selected, aes(x = Method, y = GC_Dropout, fill = Method)) +
  geom_violin(trim = FALSE, alpha = 0.8, show.legend = FALSE) +
  geom_boxplot(width=0.05, alpha = 0.6, show.legend = FALSE) +
  geom_signif(comparisons = list(c("myBaits", "Twist")),map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05),tip_length = 0,y_position=1) +
  geom_point(position = position_jitter(width = 0.2), size = 0.8, show.legend = FALSE) + 
  scale_fill_manual(values = c("myBaits" = "#E69F00", "Twist" = "#0072B2")) +
  labs(title="GC dropout",y = "", x="") +
  theme_minimal() +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8))

drop2<- ggplot(dropout_selected, aes(x = Method, y = AT_Dropout, fill = Method)) +
  geom_violin(trim = FALSE, alpha = 0.6, show.legend = FALSE) +
  geom_boxplot(width=0.05, alpha = 0.8, show.legend = FALSE) +
  geom_signif(comparisons = list(c("myBaits", "Twist")),tip_length = 0, y_position=25,map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05)) +
  geom_point(position = position_jitter(width = 0.2), size = 0.8, show.legend = FALSE) + 
  scale_fill_manual(values = c("myBaits" = "#E69F00", "Twist" = "#0072B2")) +
  labs(title="AT dropout",y = "", x="") +
  theme_minimal() +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8))

row1 <- gccov + plot_layout(guides = 'collect', widths = c(2))
row2 <- drop1 + drop2 + plot_layout(guides = 'collect', widths = c(1,1))
row1 / row2 + plot_annotation(tag_levels = 'A') & theme(plot.margin = unit(c(1, 1, 1, 1), "pt"),legend.position = "top")

# Figure 5 SNP quality ----------------------------------------------------------------------------------------------------------------
# (A) Deamination (C-T substitutions on first bp) -----------------------------------------------------------------------------------------
scale_factor <- max(dat1$X5_Prime_C.T_1st_base_on.target_104k_80bp_dedup) / max(dat1$Age,na.rm=TRUE)
deam <- ggplot(data=dat1,aes(x=Sample,y=X5_Prime_C.T_1st_base_on.target_104k_80bp_dedup, fill=Method)) + 
  geom_bar(stat="identity",position="dodge") + 
  geom_line(aes(x=Sample,y=Age* scale_factor, group = 1), color = "black", size = 0.8, linetype="dashed") +
  scale_y_continuous(name = "C-T substitutions (%)", sec.axis = sec_axis(~ . / scale_factor, name = "Sample age (ya)")) +  # Secondary axis for the rate) +
  theme_bw() + scale_fill_manual(values = c("myBaits" = "#E69F00", "Twist" = "#0072B2")) +
  scale_x_discrete(limits=mybaitsdat$Sample) +
  theme(axis.text.y.right = element_text(angle = 270, vjust = 0.5),
    axis.title.x=element_blank(),
    axis.text.x = element_text(face="bold", angle=90),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = c(0.9, 0.7),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10))

# SNP number and accuracy -------------------------------------------------------------------------------------------------
###### the filter for snpcalling on angsd is -minMapQ 20 -minQ 20 and also dp>=3 for geno calling ------------------------------
# (B) SNPs are counted as present as long as there is genotype info, which can be homozygous or heterozygous --------------
# (C) Using shotgun as the baseline, as long as one of the two alleles overlap with shotgun, the site is a match ----------
setwd("/Users/chyiyin/Dropbox/CORVID_baits/Analyses/angsd") 
mybaitsgeno <- read.table(gzfile("mybaits_104kpanel_geno_maxmis_q20_dp3.geno.gz"), sep = "\t", header=FALSE) #104K with Ns
twistgeno <- read.table(gzfile("twist_104kpanel_geno_maxmis_q20_dp3.geno.gz"), sep = "\t", header=FALSE) #104K with Ns
shotgungeno <- read.table(gzfile("shotgun_noudg_104kpanel_geno_maxmis_q20_dp3.geno.gz"), sep = "\t", header=FALSE) #104K with Ns
popshotgun_selected <- read.table("popshotgun_selected.txt", sep="\t", header=TRUE) 
twistgeno <- twistgeno[ , -ncol(twistgeno)]
mybaitsgeno <- mybaitsgeno[ , -ncol(mybaitsgeno)]

colnames(twistgeno) <-  c("chr","pos",popshotgun_selected$Samples)
colnames(mybaitsgeno) <-  c("chr","pos",popshotgun_selected$Samples)
colnames(shotgungeno) <-  c("chr","pos",popshotgun_selected$Samples)

twistgeno <- twistgeno %>% mutate(scaff_name=paste(chr,pos,sep="_")) %>% 
  column_to_rownames(var = "scaff_name") %>%
  mutate(across(everything(), ~ ifelse(. == "N", NA, .))) 
mybaitsgeno <- mybaitsgeno %>% mutate(scaff_name=paste(chr,pos,sep="_"))  %>% 
  column_to_rownames(var = "scaff_name") %>%
  mutate(across(everything(), ~ ifelse(. == "N", NA, .))) 
shotgungeno <- shotgungeno %>% mutate(scaff_name=paste(chr,pos,sep="_"))  %>% 
  column_to_rownames(var = "scaff_name") %>%
  mutate(across(everything(), ~ ifelse(. == "NN", NA, .))) 

samples <- c("BRW001", "DVT014", "NCP002", "VKP001")
lenient_match <- function(geno1, geno2) {
  if (is.na(geno1) | is.na(geno2)) return(FALSE)
  alleles1 <- strsplit(geno1, "")[[1]]
  alleles2 <- strsplit(geno2, "")[[1]]
  return(any(alleles1 %in% alleles2))
}

venn_plots <- list()
accuracy_results <- data.frame(Sample = character(), Method = character(), Accuracy = numeric())
for (sample in samples) {
  twist_sites <- rownames(twistgeno)[!is.na(twistgeno[[sample]])]
  mybaits_sites <- rownames(mybaitsgeno)[!is.na(mybaitsgeno[[sample]])]
  shotgun_sites <- rownames(shotgungeno)[!is.na(shotgungeno[[sample]])]
  
  venn_data <- list(Twist = twist_sites, myBaits = mybaits_sites, Shotgun = shotgun_sites)
  total_unique <- length(unique(unlist(venn_data)))
  venn_plot <- ggvenn(venn_data, text_size=3,show_percentage = TRUE,
    fill_color = c("#0072B2","#E69F00","darkseagreen3"), fill_alpha=0.8, stroke_size = 0.1, set_name_size = 0) +
    labs(title = paste(sample," n=",total_unique,sep="")) +
    theme_minimal()+theme_void()+
    theme(plot.title = element_text(size = 8, face = "bold", hjust = 0.5, margin = margin(b = -100)),  # ↓ Reduce bottom margin of title
          plot.margin = margin(t = -10, r = -10, b = -10, l = -10))
  venn_plots[[sample]] <- venn_plot
  
  overlap_twist <- Reduce(intersect,list(twist_sites, shotgun_sites, mybaits_sites))
  if (length(overlap_twist) > 0) {
    matches_twist <- mapply(lenient_match,
                            twistgeno[overlap_twist, sample],
                            shotgungeno[overlap_twist, sample])
    accuracy_twist <- sum(matches_twist) / length(matches_twist)
    accuracy_results <- rbind(accuracy_results, data.frame(Sample = sample, Method = "Twist", Accuracy = accuracy_twist,Overlap_Count=length(overlap_twist)))
  }
  
  overlap_mybaits <- Reduce(intersect,list(mybaits_sites, shotgun_sites, twist_sites))
  if (length(overlap_mybaits) > 0) {
    matches_mybaits <- mapply(lenient_match,
                              mybaitsgeno[overlap_mybaits, sample],
                              shotgungeno[overlap_mybaits, sample])
    accuracy_mybaits <- sum(matches_mybaits) / length(matches_mybaits)
    accuracy_results <- rbind(accuracy_results,
                              data.frame(Sample = sample, Method = "myBaits", Accuracy = accuracy_mybaits,Overlap_Count=length(overlap_mybaits)))
  }
}

accuracy_results <- accuracy_results %>%
  mutate(Accuracy_Percent = round(Accuracy * 100, 2)) %>%
  mutate(Sample = factor(Sample, levels = c("BRW001", "VKP001", "DVT014", "NCP002")))
accuracy_plot <- ggplot(accuracy_results, aes(x = Sample, y = Accuracy_Percent, fill = Method)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_text(aes(label = sprintf("%.1f", Accuracy_Percent)), position = position_dodge(width = 0.9), vjust = 1.2, size = 3) +
  geom_text(data=accuracy_results%>%filter(Method == "myBaits"), aes(x=Sample, y=0, label = Overlap_Count), vjust = 0.5, size = 3) + 
  scale_fill_manual(values = c("myBaits" = "#E69F00", "Twist" = "#0072B2","Shotgun"="darkseagreen3")) +
  labs(y = "Genotype accuracy (%)", x = "Sample") +
  theme_minimal() +theme(legend.position = "none",axis.title.x=element_blank(),panel.grid.minor = element_blank())

venn_all <- venn_plots[[samples[1]]] + venn_plots[[samples[4]]] + venn_plots[[samples[2]]] + venn_plots[[samples[3]]] + plot_layout(guides = 'collect', widths = c(1, 1,1,1))
fig5 <- deam / venn_all / accuracy_plot + plot_layout(height = c(1,1,1)) + plot_annotation(tag_levels = 'A') & theme(plot.margin = unit(c(1, 1, 1, 1), "pt")) #7x6

# Supplementary Figure 3: Allele depth for allele bias ---------------------------------------------------
##### ind0:BRW001, ind1:DVT014, ind2:NCP002, ind3:VKP001
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
  filter(sum(c(A > 0, C > 0, G > 0, T > 0)) == 2) %>% #only biallele sites
  mutate(AC = ifelse(A > 0 & C > 0, A/(A+C), 0),
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
  mutate(AC = ifelse(A > 0 & C > 0, A/(A+C), 0),
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
  mutate(AC = ifelse(A > 0 & C > 0, A/(A+C), 0),
    AG = ifelse(A > 0 & G > 0, A/(A+G), 0),
    AT = ifelse(A > 0 & T > 0, A/(A+T), 0),
    CG = ifelse(C > 0 & G > 0, C/(C+G), 0),
    CT = ifelse(C > 0 & T > 0, C/(C+T), 0),
    GT = ifelse(G > 0 & T > 0, G/(G+T), 0))

common_rows <- Reduce(intersect, list(mybaitsAD$scaff, twistAD$scaff, shotgunAD$scaff)) #11281
#common_rows <- Reduce(intersect, list(rownames(het_brw001), mybaitsAD$scaff, twistAD$scaff, shotgunAD$scaff)) #10407
mybaitsAD_filtered <- mybaitsAD[mybaitsAD$scaff %in% common_rows, ]
twistAD_filtered <- twistAD[twistAD$scaff %in% common_rows, ]
shotgunAD_filtered <- shotgunAD[shotgunAD$scaff %in% common_rows, ]

res_mybaits <- data.frame()
logistic<-function(y=0){1/(1+exp(-y))}
for (i in c("AC", "AG", "AT", "CG", "CT", "GT")) {
  allele1 <- substr(i,1,1)
  allele2 <- substr(i,2,2)
  filtered_data <- mybaitsAD_filtered %>% filter(.data[[i]] > 0)
  allele1_counts <- filtered_data[[allele1]]
  allele2_counts <- filtered_data[[allele2]]
  lm1 <- glm(cbind(allele1_counts, allele2_counts)~1,data=filtered_data,family="quasibinomial")
  p.value<-summary(lm1)$coefficients[4]
  p.value.bonf <- min(p.value * 18, 1)
  p.trans<-logistic(lm1$coefficients)
  CIlow.trans<-logistic(lm1$coefficients-summary(lm1)$coefficients[2]*1.96)
  CIup.trans<-logistic(lm1$coefficients+summary(lm1)$coefficients[2]*1.96)
  res_mybaits <- rbind(res_mybaits, data.frame(Genotype=i, estADratio=p.trans, p.value=p.value,p.value.bonf = p.value.bonf, 
                                               CI.low=CIlow.trans, CI.high=CIup.trans,type="mybaits")) }

res_twist <- data.frame()
logistic<-function(y=0){1/(1+exp(-y))}
for (i in c("AC", "AG", "AT", "CG", "CT", "GT")) {
  allele1 <- substr(i,1,1)
  allele2 <- substr(i,2,2)
  filtered_data <- twistAD_filtered %>% filter(.data[[i]] > 0)
  allele1_counts <- filtered_data[[allele1]]
  allele2_counts <- filtered_data[[allele2]]
  lm1 <- glm(cbind(allele1_counts, allele2_counts)~1,data=filtered_data,family="quasibinomial")
  p.value<-summary(lm1)$coefficients[4]
  p.value.bonf <- min(p.value * 18, 1)
  p.trans<-logistic(lm1$coefficients)
  CIlow.trans<-logistic(lm1$coefficients-summary(lm1)$coefficients[2]*1.96)
  CIup.trans<-logistic(lm1$coefficients+summary(lm1)$coefficients[2]*1.96)
  res_twist <- rbind(res_twist, data.frame(Genotype=i, estADratio=p.trans, p.value=p.value, p.value.bonf=p.value.bonf, 
                                           CI.low=CIlow.trans, CI.high=CIup.trans,type="twist")) }

res_shotgun <- data.frame()
logistic<-function(y=0){1/(1+exp(-y))}
for (i in c("AC", "AG", "AT", "CG", "CT", "GT")) {
  allele1 <- substr(i,1,1)
  allele2 <- substr(i,2,2)
  filtered_data <- shotgunAD_filtered %>% filter(.data[[i]] > 0)
  allele1_counts <- filtered_data[[allele1]]
  allele2_counts <- filtered_data[[allele2]]
  lm1 <- glm(cbind(allele1_counts, allele2_counts)~1,data=filtered_data,family="quasibinomial")
  p.value<-summary(lm1)$coefficients[4]
  p.value.bonf <- min(p.value *18, 1)
  p.trans<-logistic(lm1$coefficients)
  CIlow.trans<-logistic(lm1$coefficients-summary(lm1)$coefficients[2]*1.96)
  CIup.trans<-logistic(lm1$coefficients+summary(lm1)$coefficients[2]*1.96)
  res_shotgun <- rbind(res_shotgun, data.frame(Genotype=i, estADratio=p.trans, p.value=p.value, p.value.bonf=p.value.bonf, 
                                               CI.low=CIlow.trans, CI.high=CIup.trans,type="shotgun")) } 

abias1 <- rbind(res_shotgun, res_mybaits, res_twist) %>%
  mutate(Genotype = factor(Genotype, levels = c("AG", "CT", "AC", "AT", "CG", "GT"))) %>%
  arrange(Genotype)
  
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
  scale_fill_manual(values = c("darkseagreen3","#E69F00","#0072B2"), name = "Method", labels = c("shotgun","myBaits", "Twist")) +
  labs(y="BRW001") +
  theme_minimal() +
  theme(legend.position = "none") + 
  theme(axis.title.x = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),axis.line.x = element_blank())
a2 <- ggplot(combined_CT, aes(x = Group, y = Frequency, fill = Group)) +
  geom_violin(trim = TRUE) +  # Use trim = TRUE if you want to trim the tails
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +  
  scale_fill_manual(values = c("darkseagreen3","#E69F00","#0072B2"), name = "Method", labels = c("shotgun","myBaits", "Twist")) +
  theme_minimal() +
  theme(legend.position = "none") + 
  theme(axis.title.x = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),axis.line.x = element_blank(),axis.title.y = element_blank())
a3 <- ggplot(combined_AC, aes(x = Group, y = Frequency, fill = Group)) +
  geom_violin(trim = TRUE) +  # Use trim = TRUE if you want to trim the tails
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +  
  scale_fill_manual(values = c("darkseagreen3","#E69F00","#0072B2"), name = "Method", labels = c("shotgun","myBaits", "Twist")) +
  theme_minimal() +
  theme(legend.position = "none") + 
  theme(axis.title.x = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),axis.line.x = element_blank(),axis.title.y = element_blank())
a4 <- ggplot(combined_AT, aes(x = Group, y = Frequency, fill = Group)) +
  geom_violin(trim = TRUE) +  # Use trim = TRUE if you want to trim the tails
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +  
  scale_fill_manual(values = c("darkseagreen3","#E69F00","#0072B2"), name = "Method", labels = c("shotgun","myBaits", "Twist")) +
  theme_minimal() +
  theme(legend.position = "none") + 
  theme(axis.title.x = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),axis.line.x = element_blank(),axis.title.y = element_blank())
a5 <- ggplot(combined_CG, aes(x = Group, y = Frequency, fill = Group)) +
  geom_violin(trim = TRUE) +  # Use trim = TRUE if you want to trim the tails
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +  
  scale_fill_manual(values = c("darkseagreen3","#E69F00","#0072B2"), name = "Method", labels = c("shotgun","myBaits", "Twist")) +
  theme_minimal() +
  theme(legend.position = "none") + 
  theme(axis.title.x = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),axis.line.x = element_blank(),axis.title.y = element_blank())
a6 <- ggplot(combined_GT, aes(x = Group, y = Frequency, fill = Group)) +
  geom_violin(trim = TRUE) +  # Use trim = TRUE if you want to trim the tails
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +  
  scale_fill_manual(values = c("darkseagreen3","#E69F00","#0072B2"), name = "Method", labels = c("Shotgun","myBaits", "Twist")) +
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
    scale_fill_manual(values = c("darkseagreen3","#E69F00","#0072B2"), name = "Method", labels = c("Shotgun","myBaits", "Twist")) +
    theme_minimal() +
    theme(legend.position = "top") )

#Sample DVT014 ------------------------------------------------------------
mybaitsAD <- cbind(mybaitspos, mybaitscount[5:8]) %>%
  rename(A=4, C=5, G=6, T=7) %>%
  mutate(totDepth=A+C+G+T)  %>%
  mutate(scaff = paste(chr, pos, sep = "_")) %>%
  filter(totDepth >2) %>%
  rowwise() %>%
  filter(sum(c(A > 0, C > 0, G > 0, T > 0)) == 2) %>%
  mutate(AC = ifelse(A > 0 & C > 0, A/(A+C), 0),
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
  mutate(AC = ifelse(A > 0 & C > 0, A/(A+C), 0),
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

common_rows <- Reduce(intersect, list(mybaitsAD$scaff, twistAD$scaff, shotgunAD$scaff)) #3331
mybaitsAD_filtered <- mybaitsAD[mybaitsAD$scaff %in% common_rows, ]
twistAD_filtered <- twistAD[twistAD$scaff %in% common_rows, ]
shotgunAD_filtered <- shotgunAD[shotgunAD$scaff %in% common_rows, ]

res_mybaits <- data.frame()
logistic<-function(y=0){1/(1+exp(-y))}
for (i in c("AC", "AG", "AT", "CG", "CT", "GT")) {
  allele1 <- substr(i,1,1)
  allele2 <- substr(i,2,2)
  filtered_data <- mybaitsAD_filtered %>% filter(.data[[i]] > 0)
  allele1_counts <- filtered_data[[allele1]]
  allele2_counts <- filtered_data[[allele2]]
  lm1 <- glm(cbind(allele1_counts, allele2_counts)~1,data=filtered_data,family="quasibinomial")
  p.value<-summary(lm1)$coefficients[4]
  p.value.bonf <- min(p.value * 18, 1)
  p.trans<-logistic(lm1$coefficients)
  CIlow.trans<-logistic(lm1$coefficients-summary(lm1)$coefficients[2]*1.96)
  CIup.trans<-logistic(lm1$coefficients+summary(lm1)$coefficients[2]*1.96)
  res_mybaits <- rbind(res_mybaits, data.frame(Genotype=i, estADratio=p.trans, p.value=p.value, p.value.bonf=p.value.bonf,
                                               CI.low=CIlow.trans, CI.high=CIup.trans,type="mybaits")) }

res_twist <- data.frame()
logistic<-function(y=0){1/(1+exp(-y))}
for (i in c("AC", "AG", "AT", "CG", "CT", "GT")) {
  allele1 <- substr(i,1,1)
  allele2 <- substr(i,2,2)
  filtered_data <- twistAD_filtered %>% filter(.data[[i]] > 0)
  allele1_counts <- filtered_data[[allele1]]
  allele2_counts <- filtered_data[[allele2]]
  lm1 <- glm(cbind(allele1_counts, allele2_counts)~1,data=filtered_data,family="quasibinomial")
  p.value<-summary(lm1)$coefficients[4]
  p.value.bonf <- min(p.value * 18, 1)
  p.trans<-logistic(lm1$coefficients)
  CIlow.trans<-logistic(lm1$coefficients-summary(lm1)$coefficients[2]*1.96)
  CIup.trans<-logistic(lm1$coefficients+summary(lm1)$coefficients[2]*1.96)
  res_twist <- rbind(res_twist, data.frame(Genotype=i, estADratio=p.trans, p.value=p.value, p.value.bonf=p.value.bonf,
                                           CI.low=CIlow.trans, CI.high=CIup.trans, type="twist")) }

res_shotgun <- data.frame()
logistic<-function(y=0){1/(1+exp(-y))}
for (i in c("AC", "AG", "AT", "CG", "CT", "GT")) {
  allele1 <- substr(i,1,1)
  allele2 <- substr(i,2,2)
  filtered_data <- shotgunAD_filtered %>% filter(.data[[i]] > 0)
  allele1_counts <- filtered_data[[allele1]]
  allele2_counts <- filtered_data[[allele2]]
  lm1 <- glm(cbind(allele1_counts, allele2_counts)~1,data=filtered_data,family="quasibinomial")
  p.value<-summary(lm1)$coefficients[4]
  p.value.bonf <- min(p.value * 18, 1)
  p.trans<-logistic(lm1$coefficients)
  CIlow.trans<-logistic(lm1$coefficients-summary(lm1)$coefficients[2]*1.96)
  CIup.trans<-logistic(lm1$coefficients+summary(lm1)$coefficients[2]*1.96)
  res_shotgun <- rbind(res_shotgun, data.frame(Genotype=i, estADratio=p.trans, p.value=p.value, p.value.bonf=p.value.bonf,
                                               CI.low=CIlow.trans, CI.high=CIup.trans, type="shotgun")) }

abias2 <- rbind(res_shotgun,res_mybaits,res_twist) %>% mutate(sample = "DVT014") %>%
  mutate(Genotype = factor(Genotype, levels = c("AG", "CT", "AC", "AT", "CG", "GT"))) %>%
  arrange(Genotype)

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
  scale_fill_manual(values = c("darkseagreen3","#E69F00","#0072B2"), name = "Method", labels = c("shotgun","myBaits", "Twist")) +
  labs(y="DVT014") +
  theme_minimal() +
  theme(legend.position = "none") + 
  theme(axis.title.x = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),axis.line.x = element_blank())
b2 <- ggplot(combined_CT, aes(x = Group, y = Frequency, fill = Group)) +
  geom_violin(trim = TRUE) +  # Use trim = TRUE if you want to trim the tails
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +  
  scale_fill_manual(values = c("darkseagreen3","#E69F00","#0072B2"), name = "Method", labels = c("shotgun","myBaits", "Twist")) +
  theme_minimal() +
  theme(legend.position = "none") + 
  theme(axis.title.x = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),axis.line.x = element_blank(),axis.title.y = element_blank())
b3 <- ggplot(combined_AC, aes(x = Group, y = Frequency, fill = Group)) +
  geom_violin(trim = TRUE) +  # Use trim = TRUE if you want to trim the tails
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +  
  scale_fill_manual(values = c("darkseagreen3","#E69F00","#0072B2"), name = "Method", labels = c("shotgun","myBaits", "Twist")) +
  theme_minimal() +
  theme(legend.position = "none") + 
  theme(axis.title.x = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),axis.line.x = element_blank(),axis.title.y = element_blank())
b4 <- ggplot(combined_AT, aes(x = Group, y = Frequency, fill = Group)) +
  geom_violin(trim = TRUE) +  # Use trim = TRUE if you want to trim the tails
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +  
  scale_fill_manual(values = c("darkseagreen3","#E69F00","#0072B2"), name = "Method", labels = c("shotgun","myBaits", "Twist")) +
  theme_minimal() +
  theme(legend.position = "none") + 
  theme(axis.title.x = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),axis.line.x = element_blank(),axis.title.y = element_blank())
b5 <- ggplot(combined_CG, aes(x = Group, y = Frequency, fill = Group)) +
  geom_violin(trim = TRUE) +  # Use trim = TRUE if you want to trim the tails
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +  
  scale_fill_manual(values = c("darkseagreen3","#E69F00","#0072B2"), name = "Method", labels = c("shotgun","myBaits", "Twist")) +
  theme_minimal() +
  theme(legend.position = "none") + 
  theme(axis.title.x = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),axis.line.x = element_blank(),axis.title.y = element_blank())
b6 <- ggplot(combined_GT, aes(x = Group, y = Frequency, fill = Group)) +
  geom_violin(trim = TRUE) +  # Use trim = TRUE if you want to trim the tails
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +  
  scale_fill_manual(values = c("darkseagreen3","#E69F00","#0072B2"), name = "Method", labels = c("shotgun","myBaits", "Twist")) +
  theme_minimal() +
  theme(legend.position = "none") + 
  theme(axis.title.x = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),axis.line.x = element_blank(),axis.title.y = element_blank())

#Sample NCP002 ----------------------------------------------------------------
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
common_rows <- Reduce(intersect, list(mybaitsAD$scaff, twistAD$scaff, shotgunAD$scaff)) #72, previously 67
mybaitsAD_filtered <- mybaitsAD[mybaitsAD$scaff %in% common_rows, ]
twistAD_filtered <- twistAD[twistAD$scaff %in% common_rows, ]
shotgunAD_filtered <- shotgunAD[shotgunAD$scaff %in% common_rows, ]

res_mybaits <- data.frame()
logistic<-function(y=0){1/(1+exp(-y))}
for (i in c("AC", "AG", "AT", "CG", "CT", "GT")) {
  allele1 <- substr(i,1,1)
  allele2 <- substr(i,2,2)
  filtered_data <- mybaitsAD_filtered %>% filter(.data[[i]] > 0)
  allele1_counts <- filtered_data[[allele1]]
  allele2_counts <- filtered_data[[allele2]]
  lm1 <- glm(cbind(allele1_counts, allele2_counts)~1,data=filtered_data,family="quasibinomial")
  p.value<-summary(lm1)$coefficients[4]
  p.value.bonf <- min(p.value * 18, 1)
  p.trans<-logistic(lm1$coefficients)
  CIlow.trans<-logistic(lm1$coefficients-summary(lm1)$coefficients[2]*1.96)
  CIup.trans<-logistic(lm1$coefficients+summary(lm1)$coefficients[2]*1.96)
  res_mybaits <- rbind(res_mybaits, data.frame(Genotype=i, estADratio=p.trans, p.value=p.value,p.value.bonf=p.value.bonf,
                                               CI.low=CIlow.trans, CI.high=CIup.trans,type="mybaits")) }

res_twist <- data.frame()
logistic<-function(y=0){1/(1+exp(-y))}
for (i in c("AC", "AG", "AT", "CG", "CT", "GT")) {
  allele1 <- substr(i,1,1)
  allele2 <- substr(i,2,2)
  filtered_data <- twistAD_filtered %>% filter(.data[[i]] > 0)
  allele1_counts <- filtered_data[[allele1]]
  allele2_counts <- filtered_data[[allele2]]
  lm1 <- glm(cbind(allele1_counts, allele2_counts)~1,data=filtered_data,family="quasibinomial")
  p.value<-summary(lm1)$coefficients[4]
  p.value.bonf <- min(p.value * 18, 1)
  p.trans<-logistic(lm1$coefficients)
  CIlow.trans<-logistic(lm1$coefficients-summary(lm1)$coefficients[2]*1.96)
  CIup.trans<-logistic(lm1$coefficients+summary(lm1)$coefficients[2]*1.96)
  res_twist <- rbind(res_twist, data.frame(Genotype=i, estADratio=p.trans, p.value=p.value,p.value.bonf=p.value.bonf,
                                           CI.low=CIlow.trans, CI.high=CIup.trans,type="twist")) }

res_shotgun <- data.frame()
logistic<-function(y=0){1/(1+exp(-y))}
for (i in c("AC", "AG", "AT", "CG", "CT", "GT")) {
  allele1 <- substr(i,1,1)
  allele2 <- substr(i,2,2)
  filtered_data <- shotgunAD_filtered %>% filter(.data[[i]] > 0)
  allele1_counts <- filtered_data[[allele1]]
  allele2_counts <- filtered_data[[allele2]]
  lm1 <- glm(cbind(allele1_counts, allele2_counts)~1,data=filtered_data,family="quasibinomial")
  p.value<-summary(lm1)$coefficients[4]
  p.value.bonf <- min(p.value * 18, 1)
  p.trans<-logistic(lm1$coefficients)
  CIlow.trans<-logistic(lm1$coefficients-summary(lm1)$coefficients[2]*1.96)
  CIup.trans<-logistic(lm1$coefficients+summary(lm1)$coefficients[2]*1.96)
  res_shotgun <- rbind(res_shotgun, data.frame(Genotype=i, estADratio=p.trans, p.value=p.value,p.value.bonf=p.value.bonf,
                                               CI.low=CIlow.trans, CI.high=CIup.trans,type="shotgun")) }

abias3 <- rbind(res_shotgun,res_mybaits,res_twist) %>%mutate(sample = "NCP002") %>%
  mutate(Genotype = factor(Genotype, levels = c("AG", "CT", "AC", "AT", "CG", "GT"))) %>%
  arrange(Genotype)

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
  scale_fill_manual(values = c("darkseagreen3","#E69F00","#0072B2"), name = "Method", labels = c("shotgun","myBaits", "Twist")) +
  labs(y="NCP002") +
  theme_minimal() +
  theme(legend.position = "none") + 
  theme(axis.title.x = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),axis.line.x = element_blank())
c2 <- ggplot(combined_CT, aes(x = Group, y = Frequency, fill = Group)) +
  geom_violin(trim = TRUE) +  # Use trim = TRUE if you want to trim the tails
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +  
  scale_fill_manual(values = c("darkseagreen3","#E69F00","#0072B2"), name = "Method", labels = c("shotgun","myBaits", "Twist")) +
  theme_minimal() +
  theme(legend.position = "none") + 
  theme(axis.title.x = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),axis.line.x = element_blank(),axis.title.y = element_blank())
c3 <- ggplot(combined_AC, aes(x = Group, y = Frequency, fill = Group)) +
  geom_violin(trim = TRUE) +  # Use trim = TRUE if you want to trim the tails
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +  
  scale_fill_manual(values = c("darkseagreen3","#E69F00","#0072B2"), name = "Method", labels = c("shotgun","myBaits", "Twist")) +
  theme_minimal() +
  theme(legend.position = "none") + 
  theme(axis.title.x = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),axis.line.x = element_blank(),axis.title.y = element_blank())
c4 <- ggplot(combined_AT, aes(x = Group, y = Frequency, fill = Group)) +
  geom_violin(trim = TRUE) +  # Use trim = TRUE if you want to trim the tails
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +  
  scale_fill_manual(values = c("darkseagreen3","#E69F00","#0072B2"), name = "Method", labels = c("shotgun","myBaits", "Twist")) +
  theme_minimal() +
  theme(legend.position = "none") + 
  theme(axis.title.x = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),axis.line.x = element_blank(),axis.title.y = element_blank())
c5 <- ggplot(combined_CG, aes(x = Group, y = Frequency, fill = Group)) +
  geom_violin(trim = TRUE) +  # Use trim = TRUE if you want to trim the tails
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +  
  scale_fill_manual(values = c("darkseagreen3","#E69F00","#0072B2"), name = "Method", labels = c("shotgun","myBaits", "Twist")) +
  theme_minimal() +
  theme(legend.position = "none") + 
  theme(axis.title.x = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),axis.line.x = element_blank(),axis.title.y = element_blank())
c6 <- ggplot(combined_GT, aes(x = Group, y = Frequency, fill = Group)) +
  geom_violin(trim = TRUE) +  # Use trim = TRUE if you want to trim the tails
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +  
  scale_fill_manual(values = c("darkseagreen3","#E69F00","#0072B2"), name = "Method", labels = c("shotgun","myBaits", "Twist")) +
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
common_rows <- Reduce(intersect, list(mybaitsAD$scaff, twistAD$scaff, shotgunAD$scaff)) #9439, previous 8879
mybaitsAD_filtered <- mybaitsAD[mybaitsAD$scaff %in% common_rows, ]
twistAD_filtered <- twistAD[twistAD$scaff %in% common_rows, ]
shotgunAD_filtered <- shotgunAD[shotgunAD$scaff %in% common_rows, ]

res_mybaits <- data.frame()
logistic<-function(y=0){1/(1+exp(-y))}
for (i in c("AC", "AG", "AT", "CG", "CT", "GT")) {
  allele1 <- substr(i,1,1)
  allele2 <- substr(i,2,2)
  filtered_data <- mybaitsAD_filtered %>% filter(.data[[i]] > 0)
  allele1_counts <- filtered_data[[allele1]]
  allele2_counts <- filtered_data[[allele2]]
  lm1 <- glm(cbind(allele1_counts, allele2_counts)~1,data=filtered_data,family="quasibinomial")
  p.value<-summary(lm1)$coefficients[4]
  p.value.bonf <- min(p.value * 18, 1)
  p.trans<-logistic(lm1$coefficients)
  CIlow.trans<-logistic(lm1$coefficients-summary(lm1)$coefficients[2]*1.96)
  CIup.trans<-logistic(lm1$coefficients+summary(lm1)$coefficients[2]*1.96)
  res_mybaits <- rbind(res_mybaits, data.frame(Genotype=i, estADratio=p.trans, p.value=p.value, p.value.bonf=p.value.bonf,
                                               CI.low=CIlow.trans, CI.high=CIup.trans,type="mybaits")) }

res_twist <- data.frame()
logistic<-function(y=0){1/(1+exp(-y))}
for (i in c("AC", "AG", "AT", "CG", "CT", "GT")) {
  allele1 <- substr(i,1,1)
  allele2 <- substr(i,2,2)
  filtered_data <- twistAD_filtered %>% filter(.data[[i]] > 0)
  allele1_counts <- filtered_data[[allele1]]
  allele2_counts <- filtered_data[[allele2]]
  lm1 <- glm(cbind(allele1_counts, allele2_counts)~1,data=filtered_data,family="quasibinomial")
  p.value<-summary(lm1)$coefficients[4]
  p.value.bonf <- min(p.value * 18, 1)
  p.trans<-logistic(lm1$coefficients)
  CIlow.trans<-logistic(lm1$coefficients-summary(lm1)$coefficients[2]*1.96)
  CIup.trans<-logistic(lm1$coefficients+summary(lm1)$coefficients[2]*1.96)
  res_twist <- rbind(res_twist, data.frame(Genotype=i, estADratio=p.trans, p.value=p.value, p.value.bonf=p.value.bonf,
                                           CI.low=CIlow.trans, CI.high=CIup.trans,type="twist")) }

res_shotgun <- data.frame()
logistic<-function(y=0){1/(1+exp(-y))}
for (i in c("AC", "AG", "AT", "CG", "CT", "GT")) {
  allele1 <- substr(i,1,1)
  allele2 <- substr(i,2,2)
  filtered_data <- shotgunAD_filtered %>% filter(.data[[i]] > 0)
  allele1_counts <- filtered_data[[allele1]]
  allele2_counts <- filtered_data[[allele2]]
  lm1 <- glm(cbind(allele1_counts, allele2_counts)~1,data=filtered_data,family="quasibinomial")
  p.value<-summary(lm1)$coefficients[4]
  p.value.bonf <- min(p.value * 18, 1)
  p.trans<-logistic(lm1$coefficients)
  CIlow.trans<-logistic(lm1$coefficients-summary(lm1)$coefficients[2]*1.96)
  CIup.trans<-logistic(lm1$coefficients+summary(lm1)$coefficients[2]*1.96)
  res_shotgun <- rbind(res_shotgun, data.frame(Genotype=i, estADratio=p.trans, p.value=p.value, p.value.bonf=p.value.bonf,
                                               CI.low=CIlow.trans, CI.high=CIup.trans,type="shotgun")) }

abias4 <- rbind(res_shotgun,res_mybaits,res_twist) %>% mutate(sample = "VKP001")  %>%
  mutate(Genotype = factor(Genotype, levels = c("AG", "CT", "AC", "AT", "CG", "GT"))) %>%
  arrange(Genotype)

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
  scale_fill_manual(values = c("darkseagreen3","#E69F00","#0072B2"), name = "Method", labels = c("Shotgun","myBaits", "Twist")) +
  labs(y="VKP001") +
  theme_minimal() +
  theme(legend.position = "none") + 
  theme(axis.title.x = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),axis.line.x = element_blank())
d2 <- ggplot(combined_CT, aes(x = Group, y = Frequency, fill = Group)) +
  geom_violin(trim = TRUE) +  # Use trim = TRUE if you want to trim the tails
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +  
  scale_fill_manual(values = c("darkseagreen3","#E69F00","#0072B2"), name = "Method", labels = c("Shotgun","myBaits", "Twist")) +
  theme_minimal() +
  theme(legend.position = "none") + 
  theme(axis.title.x = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),axis.line.x = element_blank(),axis.title.y = element_blank())
d3 <- ggplot(combined_AC, aes(x = Group, y = Frequency, fill = Group)) +
  geom_violin(trim = TRUE) +  # Use trim = TRUE if you want to trim the tails
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +  
  scale_fill_manual(values = c("darkseagreen3","#E69F00","#0072B2"), name = "Method", labels = c("Shotgun","myBaits", "Twist")) +
  theme_minimal() +
  theme(legend.position = "none") + 
  theme(axis.title.x = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),axis.line.x = element_blank(),axis.title.y = element_blank())
d4 <- ggplot(combined_AT, aes(x = Group, y = Frequency, fill = Group)) +
  geom_violin(trim = TRUE) +  # Use trim = TRUE if you want to trim the tails
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +  
  scale_fill_manual(values = c("darkseagreen3","#E69F00","#0072B2"), name = "Method", labels = c("Shotgun","myBaits", "Twist")) +
  theme_minimal() +
  theme(legend.position = "none") + 
  theme(axis.title.x = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),axis.line.x = element_blank(),axis.title.y = element_blank())
d5 <- ggplot(combined_CG, aes(x = Group, y = Frequency, fill = Group)) +
  geom_violin(trim = TRUE) +  # Use trim = TRUE if you want to trim the tails
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +  
  scale_fill_manual(values = c("darkseagreen3","#E69F00","#0072B2"), name = "Method", labels = c("Shotgun","myBaits", "Twist")) +
  theme_minimal() +
  theme(legend.position = "none") + 
  theme(axis.title.x = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),axis.line.x = element_blank(),axis.title.y = element_blank())
d6 <- ggplot(combined_GT, aes(x = Group, y = Frequency, fill = Group)) +
  geom_violin(trim = TRUE) +  # Use trim = TRUE if you want to trim the tails
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +  
  scale_fill_manual(values = c("darkseagreen3","#E69F00","#0072B2"), name = "Method", labels = c("Shotgun","myBaits", "Twist")) +
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

# Supp: Performance relative to endogenous DNA ----------------------------------------------------------------------------
# Endogenous DNA from shallow shotgun seq
scatter_mybaits <- dat %>%
  select(Sample, Method, Endogenous_DNA_shotgun, Target_eff_adjusted, Target_eff_original, Ontarget_rate_adjusted, Ontarget_rate, 
         Percentage_of_mtDNA,Mean_cov_of_104k_SNPsite_incdup, Mean_cov_of_target_104k_SNPsite_dedup,Mean_cov_of_104k_80bp_incdup) %>%
  filter(Method =="myBaits")
scatter_twist <- dat %>%
  select(Sample, Method, Endogenous_DNA_shotgun, Target_eff_adjusted, Target_eff_original, Ontarget_rate_adjusted, Ontarget_rate, 
         Percentage_of_mtDNA,Mean_cov_of_104k_SNPsite_incdup, Mean_cov_of_target_104k_SNPsite_dedup,Mean_cov_of_104k_80bp_incdup) %>%
  filter(Method =="Twist")

model <- lm(Endogenous_DNA_shotgun ~ Target_eff_original, data = scatter_mybaits)
summary_model <- summary(model)
slope <- coef(model)[2]
se_slope <- summary_model$coefficients[2, 2]
t_slope <- (slope - 1) / se_slope
p_value_slope_vs_1 <- 2 * pt(-abs(t_slope), df = model$df.residual)
intercept <- coef(model)[1]
se_intercept <- summary_model$coefficients[1, 2]
t_intercept <- (intercept - 0) / se_intercept
p_value_intercept_vs_0 <- 2 * pt(-abs(t_intercept), df = model$df.residual)
p_value <- summary(model)$coefficients[2, 4]
r2 <- summary(model)$r.squared
eqn <- paste0("\nR² = ", round(r2, 3),
              "\np (slope0) = ", format.pval(p_value, digits = 3, eps = 0.001),
              "\np (slope1) = ", format.pval(p_value_slope_vs_1, digits = 3, eps = 0.001),
              "\np (int0) = ", format.pval(p_value_intercept_vs_0, digits = 3, eps = 0.001))
e3 <- ggplot(scatter_mybaits, aes(x=Endogenous_DNA_shotgun , y = Target_eff_original)) +
  geom_point(size = 1) + geom_smooth(method = "lm", se = TRUE, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  annotate("text", x = 35, y = 20, label = eqn, hjust = 0, size = 2.5) +
  labs(x = "Endogenous DNA (%)", y = "Target eff for original panel (%)", title = "myBaits") + ylim(0,100) + xlim(0,100) +
  theme_minimal() +
  theme(plot.title = element_text(size = 8, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8))

model <- lm(Endogenous_DNA_shotgun ~ Target_eff_original, data = scatter_twist)
summary_model <- summary(model)
slope <- coef(model)[2]
se_slope <- summary_model$coefficients[2, 2]
t_slope <- (slope - 1) / se_slope
p_value_slope_vs_1 <- 2 * pt(-abs(t_slope), df = model$df.residual)
intercept <- coef(model)[1]
se_intercept <- summary_model$coefficients[1, 2]
t_intercept <- (intercept - 0) / se_intercept
p_value_intercept_vs_0 <- 2 * pt(-abs(t_intercept), df = model$df.residual)
p_value <- summary(model)$coefficients[2, 4]
r2 <- summary(model)$r.squared
eqn <- paste0("\nR² = ", round(r2, 3),
              "\np (slope0) = ", format.pval(p_value, digits = 3, eps = 0.001),
              "\np (slope1) = ", format.pval(p_value_slope_vs_1, digits = 3, eps = 0.001),
              "\np (int0) = ", format.pval(p_value_intercept_vs_0, digits = 3, eps = 0.001))
e4 <- ggplot(scatter_twist, aes(x=Endogenous_DNA_shotgun , y = Target_eff_original)) +
  geom_point(size = 1) + geom_smooth(method = "lm", se = TRUE, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  annotate("text", x = 5, y = 90, label = eqn, hjust = 0, size = 2.5) +
  labs(x = "Endogenous DNA (%)", y = "Target eff for original panel (%)", title = "Twist") + ylim(0,100) + xlim(0,100) +
  theme_minimal() +
  theme(plot.title = element_text(size = 8, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8))

model <- lm(Endogenous_DNA_shotgun ~ Ontarget_rate, data = scatter_mybaits)
summary_model <- summary(model)
slope <- coef(model)[2]
se_slope <- summary_model$coefficients[2, 2]
t_slope <- (slope - 1) / se_slope
p_value_slope_vs_1 <- 2 * pt(-abs(t_slope), df = model$df.residual)
intercept <- coef(model)[1]
se_intercept <- summary_model$coefficients[1, 2]
t_intercept <- (intercept - 0) / se_intercept
p_value_intercept_vs_0 <- 2 * pt(-abs(t_intercept), df = model$df.residual)
p_value <- summary(model)$coefficients[2, 4]
r2 <- summary(model)$r.squared
eqn <- paste0("\nR² = ", round(r2, 3),
              "\np (slope0) = ", format.pval(p_value, digits = 3, eps = 0.001),
              "\np (slope1) = ", format.pval(p_value_slope_vs_1, digits = 3, eps = 0.001),
              "\np (int0) = ", format.pval(p_value_intercept_vs_0, digits = 3, eps = 0.001))
e5 <- ggplot(scatter_mybaits, aes(x=Endogenous_DNA_shotgun , y = Ontarget_rate)) +
  geom_point(size = 1) + geom_smooth(method = "lm", se = TRUE, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  annotate("text", x = 35, y = 20, label = eqn, hjust = 0, size = 2.5) +
  labs(x = "Endogenous DNA (%)", y = "On-target rate for original panel (%)", title = "myBaits") + ylim(0,100) + xlim(0,100) +
  theme_minimal() +
  theme(plot.title = element_text(size = 8, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8))

model <- lm(Endogenous_DNA_shotgun ~ Ontarget_rate, data = scatter_twist)
summary_model <- summary(model)
slope <- coef(model)[2]
se_slope <- summary_model$coefficients[2, 2]
t_slope <- (slope - 1) / se_slope
p_value_slope_vs_1 <- 2 * pt(-abs(t_slope), df = model$df.residual)
intercept <- coef(model)[1]
se_intercept <- summary_model$coefficients[1, 2]
t_intercept <- (intercept - 0) / se_intercept
p_value_intercept_vs_0 <- 2 * pt(-abs(t_intercept), df = model$df.residual)
p_value <- summary(model)$coefficients[2, 4]
r2 <- summary(model)$r.squared
eqn <- paste0("\nR² = ", round(r2, 3),
              "\np (slope0) = ", format.pval(p_value, digits = 3, eps = 0.001),
              "\np (slope1) = ", format.pval(p_value_slope_vs_1, digits = 3, eps = 0.001),
              "\np (int0) = ", format.pval(p_value_intercept_vs_0, digits = 3, eps = 0.001))
e6 <- ggplot(scatter_twist, aes(x=Endogenous_DNA_shotgun , y = Ontarget_rate)) +
  geom_point(size = 1) + geom_smooth(method = "lm", se = TRUE, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  annotate("text", x = 5, y = 90, label = eqn, hjust = 0, size = 2.5) +
  labs(x = "Endogenous DNA (%)", y = "On-target rate for original panel (%)", title = "Twist") + ylim(0,100) + xlim(0,100) +
  theme_minimal() +
  theme(plot.title = element_text(size = 8, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8))

model <- lm(Endogenous_DNA_shotgun ~ Percentage_of_mtDNA, data = scatter_mybaits)
summary_model <- summary(model)
slope <- coef(model)[2]
se_slope <- summary_model$coefficients[2, 2]
t_slope <- (slope - 1) / se_slope
p_value_slope_vs_1 <- 2 * pt(-abs(t_slope), df = model$df.residual)
intercept <- coef(model)[1]
se_intercept <- summary_model$coefficients[1, 2]
t_intercept <- (intercept - 0) / se_intercept
p_value_intercept_vs_0 <- 2 * pt(-abs(t_intercept), df = model$df.residual)
p_value <- summary(model)$coefficients[2, 4]
r2 <- summary(model)$r.squared
eqn <- paste0("\nR² = ", round(r2, 3),
              "\np (slope0) = ", format.pval(p_value, digits = 3, eps = 0.001),
              "\np (slope1) = ", format.pval(p_value_slope_vs_1, digits = 3, eps = 0.001),
              "\np (int0) = ", format.pval(p_value_intercept_vs_0, digits = 3, eps = 0.001))
e7 <- ggplot(scatter_mybaits, aes(x=Endogenous_DNA_shotgun , y = Percentage_of_mtDNA)) +
  geom_point(size = 1) + geom_smooth(method = "lm", se = TRUE, color = "black") +
  annotate("text", x = 35, y = 0.09, label = eqn, hjust = 0, size = 2.5) +
  labs(x = "Endogenous DNA (%)", y = "Mitochondrial reads (%)", title = "myBaits") + ylim(0,0.1) + xlim(0,100) +
  theme_minimal() +
  theme(plot.title = element_text(size = 8, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8))

model <- lm(Endogenous_DNA_shotgun ~ Percentage_of_mtDNA, data = scatter_twist)
summary_model <- summary(model)
slope <- coef(model)[2]
se_slope <- summary_model$coefficients[2, 2]
t_slope <- (slope - 1) / se_slope
p_value_slope_vs_1 <- 2 * pt(-abs(t_slope), df = model$df.residual)
intercept <- coef(model)[1]
se_intercept <- summary_model$coefficients[1, 2]
t_intercept <- (intercept - 0) / se_intercept
p_value_intercept_vs_0 <- 2 * pt(-abs(t_intercept), df = model$df.residual)
p_value <- summary(model)$coefficients[2, 4]
r2 <- summary(model)$r.squared
eqn <- paste0("\nR² = ", round(r2, 3),
              "\np (slope0) = ", format.pval(p_value, digits = 3, eps = 0.001),
              "\np (slope1) = ", format.pval(p_value_slope_vs_1, digits = 3, eps = 0.001),
              "\np (int0) = ", format.pval(p_value_intercept_vs_0, digits = 3, eps = 0.001))
e8 <- ggplot(scatter_twist, aes(x=Endogenous_DNA_shotgun , y = Percentage_of_mtDNA)) +
  geom_point(size = 1) + geom_smooth(method = "lm", se = TRUE, color = "black") +
  annotate("text", x = 35, y = 0.09, label = eqn, hjust = 0, size = 2.5) +
  labs(x = "Endogenous DNA (%)", y = "Mitochondrial reads (%)", title = "Twist") + ylim(0,0.1) + xlim(0,100) +
  theme_minimal() +
  theme(plot.title = element_text(size = 8, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8))

model <- lm(Endogenous_DNA_shotgun ~ Mean_cov_of_104k_80bp_incdup, data = scatter_mybaits)
summary_model <- summary(model)
slope <- coef(model)[2]
se_slope <- summary_model$coefficients[2, 2]
t_slope <- (slope - 1) / se_slope
p_value_slope_vs_1 <- 2 * pt(-abs(t_slope), df = model$df.residual)
intercept <- coef(model)[1]
se_intercept <- summary_model$coefficients[1, 2]
t_intercept <- (intercept - 0) / se_intercept
p_value_intercept_vs_0 <- 2 * pt(-abs(t_intercept), df = model$df.residual)
p_value <- summary(model)$coefficients[2, 4]
r2 <- summary(model)$r.squared
eqn <- paste0("\nR² = ", round(r2, 3),
              "\np (slope0) = ", format.pval(p_value, digits = 3, eps = 0.001),
              "\np (slope1) = ", format.pval(p_value_slope_vs_1, digits = 3, eps = 0.001),
              "\np (int0) = ", format.pval(p_value_intercept_vs_0, digits = 3, eps = 0.001))
e9 <- ggplot(scatter_mybaits, aes(x=Endogenous_DNA_shotgun , y = Mean_cov_of_104k_80bp_incdup)) +
  geom_point(size = 1) + geom_smooth(method = "lm", se = TRUE, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  annotate("text", x = 35, y = 20, label = eqn, hjust = 0, size = 2.5) +
  labs(x = "Endogenous DNA (%)", y = "Cov of comparable panel (inc dup)", title = "myBaits") + ylim(0,100) + xlim(0,100) +
  theme_minimal() +
  theme(plot.title = element_text(size = 8, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8))

model <- lm(Endogenous_DNA_shotgun ~ Mean_cov_of_104k_80bp_incdup, data = scatter_twist)
summary_model <- summary(model)
slope <- coef(model)[2]
se_slope <- summary_model$coefficients[2, 2]
t_slope <- (slope - 1) / se_slope
p_value_slope_vs_1 <- 2 * pt(-abs(t_slope), df = model$df.residual)
intercept <- coef(model)[1]
se_intercept <- summary_model$coefficients[1, 2]
t_intercept <- (intercept - 0) / se_intercept
p_value_intercept_vs_0 <- 2 * pt(-abs(t_intercept), df = model$df.residual)
p_value <- summary(model)$coefficients[2, 4]
r2 <- summary(model)$r.squared
eqn <- paste0("\nR² = ", round(r2, 3),
              "\np (slope0) = ", format.pval(p_value, digits = 3, eps = 0.001),
              "\np (slope1) = ", format.pval(p_value_slope_vs_1, digits = 3, eps = 0.001),
              "\np (int0) = ", format.pval(p_value_intercept_vs_0, digits = 3, eps = 0.001))
e10 <- ggplot(scatter_twist, aes(x=Endogenous_DNA_shotgun , y = Mean_cov_of_104k_80bp_incdup)) +
  geom_point(size = 1) + geom_smooth(method = "lm", se = TRUE, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  annotate("text", x = 5, y = 90, label = eqn, hjust = 0, size = 2.5) +
  labs(x = "Endogenous DNA (%)", y = "Cov of comparable panel (inc dup)", title = "Twist") + ylim(0,100) + xlim(0,100) +
  theme_minimal() +
  theme(plot.title = element_text(size = 8, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8))

row1e <- e3 + e4 + e5 + e6 + plot_layout(guides = 'collect', widths = c(1, 1, 1,1))
row2e <- e7 + e8 + e9 + e10 + plot_layout(guides = 'collect', widths = c(1,1,1,1))
row1e / row2e + plot_annotation(tag_levels = 'A') & theme(plot.margin = unit(c(1, 1, 1, 1), "pt"),legend.position = "bottom") #8x5

# Supp: Means and standard deviation -------------------------------------------------------------------------
dat_summary <- dat %>% group_by(Method) %>%
  summarise(across(where(is.numeric), list(Mean = ~mean(.x, na.rm = TRUE), SD = ~sd(.x, na.rm = TRUE)))) %>% 
  mutate(Method= factor(Method, levels=c("myBaits","Twist"))) %>%
  arrange(Method)

dat2 <- read_xlsx("~/Dropbox/CORVID_baits/Analyses/rawdata_allsamples/baitscomparison.xlsx") %>%
  filter(Country %in% c("PL","B") & !Sample %in% c("DVT016","DVT022","KZR002")) %>%
  mutate(All_mapped_reads_to_original_panel = ifelse(Method == "myBaits", All_mapped_reads_to_104k_121bp, All_mapped_reads_to_232k_80bp)) %>%
  mutate(Nr_comparable_mapped_reads = Nr_mapped_reads - All_non_comparable_mapped_reads) %>%
  mutate(Target_eff_adjusted = All_mapped_reads_to_104k_80bp / Nr_comparable_mapped_reads*100) %>%
  mutate(Target_eff_original = All_mapped_reads_to_original_panel / Nr_mapped_reads*100) %>%
  mutate(Target_eff_adjusted_more = ifelse(Method == "Twist", (All_mapped_reads_to_104k_80bp / (Nr_comparable_mapped_reads*104/232))*100, Target_eff_adjusted)) %>%
  mutate(Target_eff_adjusted_more2 = ifelse(Method == "Twist", (All_mapped_reads_to_104k_80bp*2 / (Nr_comparable_mapped_reads*104/232))*100, Target_eff_adjusted)) %>%
  mutate(Target_eff_original_2x = ifelse(Method == "Twist", (All_mapped_reads_to_original_panel*2 / Nr_mapped_reads*100), Target_eff_original)) %>%
  mutate(Ontarget_rate_adjusted = All_mapped_reads_to_104k_80bp / (Nr_rawreads - All_non_comparable_mapped_reads)*100) %>%
  mutate(Ontarget_rate = All_mapped_reads_to_original_panel / Nr_rawreads*100) %>%
  mutate(Obs_exp_cov = Mean_cov_of_104k_SNPsite_incdup / Expected_genomic_coverage_of_input) %>%
  mutate(Percentage_of_mtDNA = All_MT_reads / Nr_mapped_reads * 100) %>%
  mutate(Mapped_reads_to_original_panel_dedup = ifelse(Method == "myBaits", Reads_aligned_to_104k_121bp_dedup, Reads_aligned_to_232k_80bp_dedup)) %>%
  mutate(Target_eff_adjusted_dedup = Reads_aligned_to_104k_80bp_bait_dedup / (Nr_dedup_mapped_reads-Reads_to_be_ommitted_dedup)*100) %>%
  mutate(Target_eff_original_dedup = Mapped_reads_to_original_panel_dedup / Nr_dedup_mapped_reads*100) %>%
  mutate(Ontarget_rate_adjusted_dedup = Reads_aligned_to_104k_80bp_bait_dedup / (Nr_dedup_mapped_reads - Reads_to_be_ommitted_dedup)*100) %>%
  mutate(Ontarget_rate_dedup = Mapped_reads_to_original_panel_dedup / Nr_dedup_mapped_reads*100) %>%
  mutate(Percentage_of_mtDNA_dedup = MT_reads_dedup /Nr_dedup_mapped_reads * 100)

dat_summary2 <- dat2 %>% group_by(Method) %>%
  summarise(across(where(is.numeric), list(Mean = ~mean(.x, na.rm = TRUE), SD = ~sd(.x, na.rm = TRUE)))) %>% 
  mutate(Method= factor(Method, levels=c("myBaits","Twist"))) %>%
  arrange(Method)
