library(ggplot2)
library(readr)
library(dplyr)
library(plyr)

setwd("./PD_GWAS_SouthAfrica/PD_GWAS_SumStats/Step2025/outputSummaryStats")

# Conventional GWAS
# Filtered by only lead SNPs from previous studies
df1 <- read.delim("Combined_SummaryStats_LeadSNPs.tsv", header = T, sep = "\t")

# define effect directions as aligned or opposed betas
df1$direction <- ifelse(df1$beta_ext * df1$beta_Step > 0, "Aligned", "Opposite")

# filter df1 based on for SNPs that were lead SNPs in previous studies
df1_external <- df1 %>% filter(df1$leadInStep.=="no")

# classify variants as either not significant in SAPDGC or nominally significant with a p<0.05
df1_external$pval_type <- "Not Significant"
df1_external$pval_type[df1_external$p.val_Step < 0.05] <- "Nominally Significant (p<0.05)"

# plot
ggplot(df1_external, aes(beta_Step,beta_ext)) + 
  annotate('rect', xmin=0, xmax=2, ymin=-.7, ymax=0, alpha=.2, fill='pink') +
  annotate('rect', xmin=-1.5, xmax=0, ymin=-.7, ymax=0, alpha=.2, fill='palegreen') +
  annotate('rect', xmin=-1.5, xmax=0, ymin=0, ymax=.75, alpha=.2, fill='pink') +
  annotate('rect', xmin=0, xmax=2, ymin=0, ymax=.75, alpha=.2, fill='palegreen') +
  geom_point(aes(color=pval_type)) +
  scale_color_manual(values = c( "Not Significant" = "grey60", "Nominally Significant (p<0.05)" = "deepskyblue")) +
  labs(title = "Beta-Beta Plot of Previously Identified Lead SNPs",
       x = "Beta Values for SAPDSC",
       y = "Beta Values for Previous GWAS") +
  theme(legend.title = element_blank()) +
  facet_wrap(~ study)
ggsave("./plots/PrevSigPeakVars_bySignificance.png",width=16, height = 6, dpi=600)

# Filtered by only lead SNPs from SAPDGC
# filter df1 based on for SNPs that were lead SNPs in previous studies
filteredInternaldf1 <- df1 %>% filter(df1$leadInStep. == "yes")

# filter based on variants with p<0.05 in Loesch et al. 2021
internal_df1 <- filteredInternaldf1 %>% filter(p.val_ext < 0.05)

# filtering for SNPs from Loesch et al. 2021 with p-val < 0.05
loesch_df1 <- filteredInternaldf1 %>% filter(filteredInternaldf1$study == "Loesch et al. 2021")
filtered_loesch_df1 <- loesch_df1 %>% filter(loesch_df1$p.val_ext < 0.05)

# plot
ggplot(filtered_loesch_df1, aes(beta_Step,beta_ext,color=study)) + 
  annotate('rect', xmin=0, xmax=3, ymin=-1.5, ymax=0, alpha=.2, fill='pink') +
  annotate('rect', xmin=-3, xmax=0, ymin=-1.5, ymax=0, alpha=.2, fill='palegreen') +
  annotate('rect', xmin=-3, xmax=0, ymin=0, ymax=1.5, alpha=.2, fill='pink') +
  annotate('rect', xmin=0, xmax=3, ymin=0, ymax=1.5, alpha=.2, fill='palegreen') +
  geom_point(color = "black", alpha = 0.5) +
  labs(title = "Beta-Beta Plot of Suggestively Significant lead SNPs in SAPDSC (p<0.05 in LARGE-PD)",
       x = "Beta Values for SAPDSC",
       y = "Beta Values for LARGE-PD") +
  theme(legend.title = element_blank(), plot.title = element_text(size = 10))
ggsave("./plots/StepSuggestiveSigPeakVarsAndLPDSig.png",dpi=600)

# LA-GWAS with ATT
# Replication of lead SNPs with Loesch et al. 2021
df2 <- read.delim("ATT_Combined_SummaryStats.tsv", header = T, sep = "\t")
# define effect directions as aligned or opposed betas
df2$direction <- ifelse(df2$beta_ext * df2$beta_Step > 0, "Aligned", "Opposite")

plotTitle = "ANC Beta-Beta Plot of Previously Identified Lead SNPs (p<0.05 in SAPDSC)"

# filter by external p-val < 0.05
df2_filtered <- df2 %>% filter(df2$p.val_ext<0.05)

df2_loesch <- df2_filtered %>% filter(study == "Loesch et al. 2021")
df2_loesch <- df2_loesch %>% filter(df2_loesch$leadInStep.=="yes")

# plot
ggplot(df2_loesch, aes(beta_Step,beta_ext)) + 
  annotate('rect', xmin=0, xmax=1.5, ymin=-.7, ymax=0, alpha=.2, fill='pink') +
  annotate('rect', xmin=-1.5, xmax=0, ymin=-.7, ymax=0, alpha=.2, fill='palegreen') +
  annotate('rect', xmin=-1.5, xmax=0, ymin=0, ymax=.7, alpha=.2, fill='pink') +
  annotate('rect', xmin=0, xmax=1.5, ymin=0, ymax=.7, alpha=.2, fill='palegreen') +
  geom_point(color = "black", alpha = 0.5) +
  labs(title = "Beta-Beta Plot of Suggestively Significant Variants (p<0.05 in LARGE-PD) for LA-GWAS with ATT",
       x = "Beta Values for SAPDSC",
       y = "Beta Values for LARGE-PD") +
  theme(legend.title = element_blank(), plot.title = element_text(size = 10))
ggsave("./plots/StepSuggestive_LPD_ATT.png",width=9, height = 6, dpi=600)

# filter df1 based on for SNPs that were lead SNPs in previous studies
df2_external <- df2 %>% filter(df2$leadInStep.=="no")
df2_external <- df2_external %>% filter(df2_external$study != "Loesch et al. 2021")

# classify variants as either not significant in SAPDGC or nominally significant with a p<0.05
df2_external$pval_type <- "Not Significant"
df2_external$pval_type[df2_external$p.val_Step < 0.05] <- "Nominally Significant (p<0.05)"

