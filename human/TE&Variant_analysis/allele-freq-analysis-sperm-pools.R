# =========================
# Load data
# =========================
df <- read.table("SU_sig.tsv", header=TRUE, sep="\t", stringsAsFactors=FALSE)

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
cellage <- scan("cellage_list.txt", what="character", quiet=TRUE)
genage  <- scan("genage_list.txt",  what="character", quiet=TRUE)
senmayo <- scan("senmayo_list.txt", what="character", quiet=TRUE)

# Assign category
df$Category <- "Other"
df$Category[df$Gene %in% cellage] <- "CellAge"
df$Category[df$Gene %in% genage]  <- "GenAge"
df$Category[df$Gene %in% senmayo] <- "SenMayo"

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
run_mwu("GenAge")
run_mwu("SenMayo")


# =========================
# 4. PLOTS
# =========================

# Boxplot by donor
png("delta_by_donor.png", width=800, height=600)
boxplot(Delta ~ Donor, data=df, main="Delta by Donor")
abline(h=0, lty=2)
dev.off()

# Boxplot by category
png("delta_by_category.png", width=800, height=600)
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
df <- read.table("SU_with_delta.tsv", header=TRUE, sep="\t", stringsAsFactors=FALSE)

# =========================
# Filter SIGNIFICANT variants
# =========================
df <- df[df$LRT > 10, ]

cat("Number of significant variants:", nrow(df), "\n\n")

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
write.table(donor_summary, "donor_results.tsv", sep="\t", quote=FALSE, row.names=FALSE)
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

write.table(dc_summary, "donor_comparison_results.tsv",
            sep="\t", quote=FALSE, row.names=FALSE)


# =========================
# 4. PLOT (publication-quality)
# =========================
png("delta_by_donor.png", width=900, height=700)

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
write.table(results, "gene_set_donor_comparison.tsv",
            sep="\t", quote=FALSE, row.names=FALSE)

print(results)
