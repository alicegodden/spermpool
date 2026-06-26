# =========================
# Load data
# =========================

setwd("~/Desktop/Projects/ageing/age/new_13_invesgtigation")

#DM_SU_with_delta.tsv made using DM human_sperm_sig_dec23 spreadsheet- pasted into terminal, and then
#awk 'BEGIN{FS=OFS="\t"} NR==1{print $0,"Delta"; next} {print $0,$9-$8}' SU.tsv > SU_with_delta.tsv
df <- read.table("DM_SU_with_delta.tsv", header=TRUE, sep="\t", stringsAsFactors=FALSE)

# -------------------------
# Extract donor + comparison
# -------------------------
df$Donor <- sub("\\(.*", "", df$Sample)
df$Comp  <- sub(".*\\(|\\)", "", df$Sample)

# =========================
# 1. GLOBAL TESTS
# =========================

cat("=== GLOBAL TESTS ===\n")

# Sign test
neg <- sum(df$Delta < 0)
n   <- nrow(df)
cat("Proportion negative:", neg/n, "\n")

sign_test <- binom.test(neg, n, p=0.5, alternative="greater")
print(sign_test)

# Wilcoxon signed-rank test
wilcox_global <- wilcox.test(df$Delta, mu=0, alternative="less")
print(wilcox_global)


# =========================
# 2. PER DONOR × COMPARISON
# =========================
library(dplyr)

cat("\n=== PER DONOR × COMPARISON ===\n")

summary_stats <- df %>%
  group_by(Donor, Comp) %>%
  summarise(
    n = n(),
    neg = sum(Delta < 0),
    prop = neg / n,
    .groups = "drop"
  )

# Add binomial p-values
summary_stats$p_sign <- mapply(function(neg, n) {
  binom.test(neg, n, p=0.5, alternative="greater")$p.value
}, summary_stats$neg, summary_stats$n)

# Add Wilcoxon p-values
wilcox_p <- df %>%
  group_by(Donor, Comp) %>%
  summarise(
    p_wilcox = wilcox.test(Delta, mu=0, alternative="less")$p.value,
    .groups = "drop"
  )

summary_stats <- merge(summary_stats, wilcox_p, by=c("Donor","Comp"))

print(summary_stats)


# =========================
# 3. GENE SET STRATIFICATION
# =========================

cat("\n=== GENE SET ANALYSIS ===\n")

# Load gene lists
cellage <- scan("genes13_list_ids.txt", what="character", quiet=TRUE)


# Assign category
df$Category <- "Other"
df$Category[df$Gene %in% cellage] <- "CellAge"

# Function for MWU test
run_mwu <- function(cat_name) {
  g1 <- df$Delta[df$Category == cat_name]
  g2 <- df$Delta[df$Category != cat_name]
  
  test <- wilcox.test(g1, g2)
  
  cat("\n", cat_name, "vs Other:\n")
  print(test)
  
  cat("Median", cat_name, ":", median(g1), "\n")
  cat("Median Other:", median(g2), "\n")
}

run_mwu("CellAge")


# =========================
# 4. PLOTS
# =========================

# Boxplot by donor
png("delta_by_donor_13.png", width=800, height=600)
boxplot(Delta ~ Donor, data=df, main="Delta by Donor")
abline(h=0, lty=2)
dev.off()

# Boxplot by category
png("delta_by_category_13.png", width=800, height=600)
boxplot(Delta ~ Category, data=df, main="Delta by Gene Category")
abline(h=0, lty=2)
dev.off()

cat("\nPlots saved:\n")
cat(" - delta_by_donor.png\n")
cat(" - delta_by_category.png\n")

# now to generate results for the paper

# =========================
# Load libraries
# =========================
library(dplyr)

# =========================
# Load data
# =========================
#df <- read.table("SU_with_delta.tsv", header=TRUE, sep="\t", stringsAsFactors=FALSE)

# =========================
# Filter SIGNIFICANT variants
# =========================
#df <- df[df$LRT > 10, ]

#cat("Number of significant variants:", nrow(df), "\n\n")

# =========================
# Extract donor + comparison
# =========================
df$Donor <- sub("\\(.*", "", df$Sample)
df$Comp  <- sub(".*\\(|\\)", "", df$Sample)

# =========================
# 1. GLOBAL RESULTS
# =========================
cat("=== GLOBAL RESULTS ===\n")

neg <- sum(df$Delta < 0)
n   <- nrow(df)
prop <- neg / n

cat("Negative variants:", neg, "/", n, "\n")
cat("Proportion negative:", round(prop, 4), "\n\n")

# Sign test
sign_test <- binom.test(neg, n, p=0.5, alternative="greater")
cat("Sign test p-value:", sign_test$p.value, "\n")

# Wilcoxon signed-rank test
wilcox_global <- wilcox.test(df$Delta, mu=0, alternative="less")
cat("Wilcoxon p-value:", wilcox_global$p.value, "\n")

cat("Median Delta:", median(df$Delta), "\n\n")


# =========================
# 2. DONOR-LEVEL RESULTS
# =========================
cat("=== DONOR RESULTS ===\n")

donor_summary <- df %>%
  group_by(Donor) %>%
  summarise(
    n = n(),
    neg = sum(Delta < 0),
    prop = neg / n,
    .groups = "drop"
  )

# Add p-values
donor_summary$p_sign <- mapply(function(neg, n) {
  binom.test(neg, n, p=0.5, alternative="greater")$p.value
}, donor_summary$neg, donor_summary$n)

# Wilcoxon per donor
wilcox_donor <- df %>%
  group_by(Donor) %>%
  summarise(
    p_wilcox = wilcox.test(Delta, mu=0, alternative="less")$p.value,
    .groups = "drop"
  )

donor_summary <- merge(donor_summary, wilcox_donor, by="Donor")

print(donor_summary)

# Save table for supplement
write.table(donor_summary, "donor_results_13.tsv", sep="\t", quote=FALSE, row.names=FALSE)
# addition of time point analysis too


# =========================
# 3. DONOR × COMPARISON (optional but strong)
# =========================
cat("\n=== DONOR × COMPARISON ===\n")

dc_summary <- df %>%
  group_by(Donor, Comp) %>%
  summarise(
    n = n(),
    neg = sum(Delta < 0),
    prop = neg / n,
    .groups = "drop"
  )

dc_summary$p_sign <- mapply(function(neg, n) {
  binom.test(neg, n, p=0.5, alternative="greater")$p.value
}, dc_summary$neg, dc_summary$n)

print(dc_summary)

write.table(dc_summary, "donor_comparison_results_13.tsv",
            sep="\t", quote=FALSE, row.names=FALSE)


# =========================
# 4. PLOT (publication-quality)
# =========================
png("delta_by_donor_13.png", width=900, height=700)

boxplot(Delta ~ Donor, data=df,
        col="lightgrey",
        border="black",
        ylab="Delta (Maf_O - Maf_C)",
        xlab="Donor",
        main="Distribution of allele frequency change by donor")

abline(h=0, lty=2)

dev.off()

cat("\nPlot saved: delta_by_donor.png\n")
# plot time points too
# facet

library(stringr)

df$Donor <- str_extract(df$Sample, "^Donor[0-9]+")
df$Comp  <- str_extract(df$Sample, "T0xT[0-9]+")
table(df$Comp)

ggplot(df, aes(x=Donor, y=Delta, fill=Comp)) +
  geom_boxplot(position=position_dodge(width=0.8), outlier.size=0.5) +
  geom_hline(yintercept=0, linetype="dashed") +
  labs(
    title="Delta by donor across timepoints",
    x="Donor",
    y="Delta (Maf_O - Maf_C)",
    fill="Comparison"
  ) +
  theme_bw()

# now to order the timepoints
df$Comp <- factor(
  df$Comp,
  levels = c("T0xT2", "T0xT4", "T0xT8", "T0xT12", "T0xT24", "T0xT48")
)

ggplot(df, aes(x=Donor, y=Delta, fill=Comp)) +
  geom_boxplot(position=position_dodge(width=0.8), outlier.size=0.5) +
  geom_hline(yintercept=0, linetype="dashed") +
  labs(
    title="Delta by donor across timepoints",
    x="Donor",
    y="Delta (Maf_O - Maf_C)",
    fill="Comparison"
  ) +
  theme_classic()

# make it pretty


ggplot(df, aes(x=Donor, y=Delta, fill=Comp)) +
  geom_boxplot(position=position_dodge(width=0.8), outlier.size=0.5) +
  geom_hline(yintercept=0, linetype="dashed", colour="black") +
  labs(
    title=NULL,
    x="Donor",
    y="Δ allele frequency (Maf_O - Maf_C)",
    fill="Comparison"
  ) +
  scale_fill_manual(
    values = c(
      "T0xT2"  = "#deebf7",
      "T0xT4"  = "#c6dbef",
      "T0xT8"  = "#9ecae1",
      "T0xT12" = "#6baed6",
      "T0xT24" = "#3182bd",
      "T0xT48" = "#08519c"
    )
  ) +
  theme_classic() +
  theme(
    text = element_text(face="bold"),
    axis.text.x = element_text(angle=45, hjust=1),
    axis.line = element_line(colour="black"),
    panel.grid = element_blank()
  )


library(dplyr)

# =========================
# Load data
# =========================
df <- read.table("SU_with_delta.tsv", header=TRUE, sep="\t", stringsAsFactors=FALSE)

# Filter significant variants
df <- df[df$LRT > 10, ]

# Extract donor + comparison
df$Donor <- sub("\\(.*", "", df$Sample)
df$Comp  <- sub(".*\\(|\\)", "", df$Sample)

# =========================
# Load gene lists
# =========================
cellage <- scan("cellage_list.txt", what="character", quiet=TRUE)
genage  <- scan("genage_list.txt",  what="character", quiet=TRUE)
senmayo <- scan("senmayo_list.txt", what="character", quiet=TRUE)

# =========================
# Overlap-aware flags
# =========================
df$CellAge <- df$Gene %in% cellage
df$GenAge  <- df$Gene %in% genage
df$SenMayo <- df$Gene %in% senmayo

# =========================
# Function for one test
# =========================
run_test <- function(data, flag) {
  g1 <- data$Delta[flag]
  g2 <- data$Delta[!flag]
  
  # avoid errors if too few values
  if(length(g1) < 5 | length(g2) < 5) return(NA)
  
  return(wilcox.test(g1, g2)$p.value)
}

# =========================
# Apply per Donor × Comparison
# =========================
results <- df %>%
  group_by(Donor, Comp) %>%
  summarise(
    n = n(),
    
    # counts in each set
    n_cellage = sum(CellAge),
    n_genage  = sum(GenAge),
    n_senmayo = sum(SenMayo),
    
    # medians
    med_cellage = median(Delta[CellAge]),
    med_genage  = median(Delta[GenAge]),
    med_senmayo = median(Delta[SenMayo]),
    med_other   = median(Delta[!CellAge & !GenAge & !SenMayo]),
    
    # p-values
    p_cellage = run_test(cur_data(), CellAge),
    p_genage  = run_test(cur_data(), GenAge),
    p_senmayo = run_test(cur_data(), SenMayo),
    
    .groups = "drop"
  )

# Save results
write.table(results, "gene_set_donor_comparison_nolrt.tsv",
            sep="\t", quote=FALSE, row.names=FALSE)

print(results)

# now to look at the cell age stress-induced cell sensescence genes
#taking all induces and stres-induced genes for cell sensescence cellage
stress_cellage <- scan("genes13_list_ids.txt", what="character", quiet=TRUE)

df$StressCellAge <- df$Gene %in% stress_cellage


run_mwu_stress <- function() {
  g1 <- df$Delta[df$StressCellAge]
  g2 <- df$Delta[!df$StressCellAge]
  
  test <- wilcox.test(g1, g2)
  
  cat("\nStressCellAge vs Other:\n")
  print(test)
  
  cat("Median StressCellAge:", median(g1), "\n")
  cat("Median Other:", median(g2), "\n")
}

run_mwu_stress()
df$StressCellAge   # TRUE / FALSE


png("delta_stress_vs_other_nolrt.png", width=800, height=600)

boxplot(Delta ~ StressCellAge, data=df,
        names = c("Other", "Stress CellAge"),
        col   = c("grey70", "firebrick"),
        ylab  = "Δ allele frequency (Maf_O - Maf_C)",
        main  = "Stress-induced CellAge genes vs background")

abline(h = 0, lty = 2)

dev.off()


# now per donor
stress_cellage <- scan("genes13_list_ids.txt", what="character", quiet=TRUE)
df$StressCellAge <- df$Gene %in% stress_cellage


library(ggplot2)

ggplot(df, aes(x = Donor, y = Delta, fill = StressCellAge)) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~Comp) +
  scale_fill_manual(
    values = c("grey70", "firebrick"),
    labels = c("Other", "Stress CellAge")
  ) +
  labs(
    x = "Donor",
    y = "Δ allele frequency (Maf_O - Maf_C)",
    fill = "Gene set"
  ) +
  theme_classic() +
  theme(
    text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# plot in time point order
# Order timepoints properly
df$Comp <- factor(
  df$Comp,
  levels = c("T0xT2", "T0xT4", "T0xT8", "T0xT12", "T0xT24", "T0xT48")
)

# Clean any stray characters (like brackets)
df$Comp <- gsub("[()]", "", df$Comp)



library(ggplot2)

ggplot(df, aes(x = Donor, y = Delta, fill = StressCellAge)) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  
  facet_wrap(~Comp) +
  
  scale_fill_manual(
    values = c("grey70", "firebrick"),
    labels = c("Other", "Stress CellAge")
  ) +
  
  labs(
    x = "Donor",
    y = "Δ allele frequency (Maf_O - Maf_C)",
    fill = "Gene set"
  ) +
  
  theme_classic() +
  theme(
    text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold")   # <- cleaner facet titles
  )

facet_wrap(~Comp, labeller = label_value)

df$Comp <- recode(df$Comp,
                  "T0xT2"  = "T0–T2",
                  "T0xT4"  = "T0–T4",
                  "T0xT8"  = "T0–T8",
                  "T0xT12" = "T0–T12",
                  "T0xT24" = "T0–T24",
                  "T0xT48" = "T0–T48")



library(stringr)

df$Donor <- str_extract(df$Sample, "^Donor[0-9]+")
df$Comp  <- str_extract(df$Sample, "T0xT[0-9]+")

df$Comp <- factor(
  df$Comp,
  levels = c("T0xT2", "T0xT4", "T0xT8", "T0xT12", "T0xT24", "T0xT48")
)


table(df$Comp)


ggplot(df, aes(x = Donor, y = Delta, fill = StressCellAge)) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  
  facet_wrap(~Comp, nrow = 1) +
  
  scale_fill_manual(
    values = c("grey70", "firebrick"),
    labels = c("Other", "Stress CellAge")
  ) +
  
  labs(
    x = "Donor",
    y = "Δ allele frequency (Maf_O - Maf_C)"
  ) +
  
  theme_classic() +
  theme(
    # ✅ Make everything bold
    text = element_text(face = "bold"),
    
    # ✅ Rotate donor labels
    axis.text.x = element_text(angle = 45, hjust = 1),
    
    # ✅ Keep facet labels bold and clear
    strip.text = element_text(face = "bold"),
    
    # ✅ Clean styling
    axis.line = element_line(colour = "black"),
    panel.grid = element_blank()
  )

#stats

wilcox.test(df$Delta[df$StressCellAge],
            df$Delta[!df$StressCellAge])

median(df$Delta[df$StressCellAge])
median(df$Delta[!df$StressCellAge])


# filter to the list of 12
#gene	ensembl_id
#AKR1B1	ENSG00000085662
#ALDH2	ENSG00000111275
#ATF7IP	ENSG00000171681
#HAS1	ENSG00000105723
#HS2ST1	ENSG00000171476
#KCNJ12	ENSG00000171223
#ME1	ENSG00000014641
#NBN	ENSG00000160877
#NTN4	ENSG00000213281
#PRKDC	ENSG00000162409
#TACC3	ENSG00000188612
#ZDHHC3	ENSG00000137693


ids_of_interest <- c(
  "ENSG00000085662", "ENSG00000111275", "ENSG00000171681",
  "ENSG00000105723", "ENSG00000171476", "ENSG00000171223",
  "ENSG00000014641", "ENSG00000160877", "ENSG00000213281",
  "ENSG00000162409", "ENSG00000188612", "ENSG00000137693"
)

df_filtered <- df %>%
  filter(
    Gene %in% ids_of_interest,
    StressCellAge == "Stress CellAge"
  )


ggplot(df_filtered, aes(x = Donor, y = Delta, fill = StressCellAge)) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~Comp, nrow = 1) +
  scale_fill_manual(
    values = c("grey70", "firebrick"),
    labels = c("Other", "Stress CellAge")
  ) +
  labs(
    x = "Donor",
    y = "Δ allele frequency (Maf_O - Maf_C)"
  ) +
  theme_classic() +
  theme(
    text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold"),
    axis.line = element_line(colour = "black"),
    panel.grid = element_blank()
  )


df <- df %>%
  mutate(Gene_clean = sub("\\..*", "", Gene))
df_filtered <- df %>%
  filter(
    Gene_clean %in% ids_of_interest,
    StressCellAge == "Stress CellAge"
  )

df %>% filter(Gene_clean %in% ids_of_interest) %>% nrow()

df <- df %>%
  mutate(
    Gene_clean = sub("\\..*", "", Gene),
    GeneSet = ifelse(Gene_clean %in% ids_of_interest,
                     "Selected genes", "Other")
  )

ggplot(df, aes(x = Donor, y = Delta, fill = GeneSet)) +
  geom_boxplot(outlier.size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~Comp, nrow = 1) +
  scale_fill_manual(values = c("grey80", "firebrick")) +
  theme_classic() +
  theme(
    text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


# to filter to plot just the donors 5-9
library(dplyr)

df %>%
  filter(Donor %in% paste0("Donor", 5:9)) %>%
  ggplot(aes(x = Donor, y = Delta, fill = StressCellAge)) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~Comp) +
  scale_fill_manual(
    values = c("grey70", "firebrick"),
    labels = c("Other", "Stress CellAge")
  ) +
  labs(
    x = "Donor",
    y = "Δ allele frequency (Maf_O - Maf_C)",
    fill = "Gene set"
  ) +
  theme_classic() +
  theme(
    text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold")
  )

# now just a few timepoints
library(dplyr)

df %>%
  filter(
    Donor %in% paste0("Donor", 5:9),
    tolower(Comp) %in% c("t0xt4", "t0xt24", "t0xt48")
  ) %>%
  mutate(
    Donor = factor(Donor, levels = paste0("Donor", 5:9))
  ) %>%
  ggplot(aes(x = Donor, y = Delta, fill = StressCellAge)) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~Comp) +
  scale_fill_manual(
    values = c("grey70", "firebrick"),
    labels = c("Other", "Stress CellAge")
  ) +
  labs(
    x = "Donor",
    y = "Δ allele frequency (Maf_O - Maf_C)",
    fill = "Gene set"
  ) +
  theme_classic() +
  theme(
    text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold")
  )


# how many of the "13" made it
genes13 <- unique(read.table("genes13_list_ids.txt", stringsAsFactors = FALSE)[,1])
sum(genes13 %in% df$Gene)
intersect(genes13, df$Gene)
