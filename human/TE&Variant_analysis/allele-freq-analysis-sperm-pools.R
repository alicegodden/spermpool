# ==============================================================
# Title:  Analysis of Delta Allele Frequency Changes
# Author: Dr. Alice M. Godden
# ==============================================================

# =========================
# 0. Setup
# =========================

library(dplyr)
library(ggplot2)
library(stringr)
library(data.table)

# Set working directory (EDIT PATH AS NEEDED)
setwd("~/Desktop/Projects/ageing/age/new_13_invesgtigation")

# =========================
# 1. Helper Functions
# =========================

# Extract donor + comparison from Sample column
extract_metadata <- function(df) {
  df %>%
    mutate(
      Donor = sub("\\(.*", "", Sample),
      Comp  = sub(".*\\(|\\)", "", Sample)
    )
}

# Run sign test
run_sign_test <- function(neg, n) {
  binom.test(neg, n, p = 0.5, alternative = "greater")$p.value
}

# Safe Wilcoxon test (avoids crashes)
safe_wilcox <- function(x, y = NULL, mu = 0, paired = FALSE, alt = "less") {
  tryCatch({
    if (is.null(y)) {
      wilcox.test(x, mu = mu, alternative = alt)$p.value
    } else {
      if (length(x) < 5 || length(y) < 5) return(NA)
      wilcox.test(x, y)$p.value
    }
  }, error = function(e) NA)
}

# Generic grouped summary
summarise_delta <- function(df, group_vars) {
  df %>%
    group_by(across(all_of(group_vars))) %>%
    summarise(
      n = n(),
      neg = sum(Delta < 0),
      prop = neg / n,
      p_sign = run_sign_test(neg, n),
      p_wilcox = safe_wilcox(Delta),
      .groups = "drop"
    )
}

# MWU test for gene sets
run_mwu <- function(df, flag, label) {
  g1 <- df$Delta[flag]
  g2 <- df$Delta[!flag]

  cat("\n---", label, "vs Other ---\n")
  cat("Median", label, ":", median(g1, na.rm = TRUE), "\n")
  cat("Median Other:", median(g2, na.rm = TRUE), "\n")
  cat("p-value:", safe_wilcox(g1, g2), "\n")
}

# =========================
# 2. Load Data
# =========================

df <- fread("DM_SU_with_delta.tsv") %>%
  as.data.frame()

# Optional filter (toggle as needed)
# df <- df[df$LRT > 10, ]

# Extract metadata
df <- extract_metadata(df)

# =========================
# 3. Global Tests
# =========================

cat("=== GLOBAL TESTS ===\n")

neg <- sum(df$Delta < 0)
n   <- nrow(df)

cat("Proportion negative:", neg/n, "\n")
cat("Sign test p:", run_sign_test(neg, n), "\n")
cat("Wilcoxon p:", safe_wilcox(df$Delta), "\n")
cat("Median Delta:", median(df$Delta), "\n")

# =========================
# 4. Grouped Results
# =========================

cat("\n=== DONOR RESULTS ===\n")
donor_summary <- summarise_delta(df, "Donor")
print(donor_summary)

write.table(donor_summary, "donor_results.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("\n=== DONOR × COMPARISON ===\n")
dc_summary <- summarise_delta(df, c("Donor", "Comp"))
print(dc_summary)

write.table(dc_summary, "donor_comparison_results.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

# =========================
# 5. Gene Set Analysis
# =========================

# Load gene sets
cellage <- scan("genes13_list_ids.txt", what = "character", quiet = TRUE)

# Annotate
df <- df %>%
  mutate(CellAge = Gene %in% cellage)

run_mwu(df, df$CellAge, "CellAge")

# =========================
# 6. Plotting Helpers
# =========================

plot_box <- function(data, formula, file, title, ylab = "Delta") {
  png(file, width = 800, height = 600)
  boxplot(formula,
          data = data,
          main = title,
          ylab = ylab,
          col = "lightgrey",
          border = "black")
  abline(h = 0, lty = 2)
  dev.off()
}

# Basic plots
plot_box(df, Delta ~ Donor, "delta_by_donor.png", "Delta by Donor")
plot_box(df, Delta ~ CellAge, "delta_by_category.png", "Delta by Gene Category")

# =========================
# 7. Timepoint Plot (ggplot)
# =========================

# Clean + order timepoints
df <- df %>%
  mutate(
    Comp = str_extract(Sample, "T0xT[0-9]+"),
    Comp = factor(Comp,
                  levels = c("T0xT2","T0xT4","T0xT8","T0xT12","T0xT24","T0xT48"))
  )

p_time <- ggplot(df, aes(x = Donor, y = Delta, fill = Comp)) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_classic() +
  labs(x = "Donor",
       y = "Δ allele frequency",
       fill = "Comparison") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save

ggsave("delta_timepoint_plot.png", p_time, width = 10, height = 6, dpi = 300)

# =========================
# 8. VEP + CADD Integration
# =========================

parse_csq <- function(info) {
  csq <- str_extract(info, "CSQ=.*")
  csq <- sub("CSQ=", "", csq)

  first <- str_split(csq, ",")[[1]][1]
  fields <- str_split(first, "\\|")[[1]]

  nums <- suppressWarnings(as.numeric(fields))
  nums <- nums[!is.na(nums)]

  data.frame(
    IMPACT = fields[3],
    CADD = ifelse(length(nums) > 0, max(nums), NA)
  )
}

# Load VEP files
vcf <- fread("vep_input.txt", skip = "#CHROM")
colnames(vcf) <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO")

parsed <- do.call(rbind, lapply(vcf$INFO, parse_csq))
vcf <- cbind(vcf, parsed)

# Merge key
make_key <- function(chr, pos, ref, alt) {
  paste(chr, pos, ref, alt, sep = "_")
}

df$key  <- make_key(df$Chromosome, df$ranges, df$Ref, df$Alt)
vcf$key <- make_key(vcf$CHROM, vcf$POS, vcf$REF, vcf$ALT)

merged <- left_join(df, vcf[, c("key","CADD","IMPACT")], by = "key")

# =========================
# 9. CADD Plot
# =========================

plot_df <- merged %>%
  filter(!is.na(CADD)) %>%
  mutate(
    CADD_cat = case_when(
      CADD < 10 ~ "Benign",
      CADD < 20 ~ "Moderate",
      CADD < 30 ~ "Deleterious",
      TRUE ~ "Highly deleterious"
    ),
    direction = ifelse(Delta > 0, "Increase", "Decrease")
  )

p_cadd <- ggplot(plot_df, aes(x = Gene, y = Delta)) +
  geom_point(aes(size = CADD,
                 colour = CADD_cat,
                 shape = direction),
             alpha = 0.8,
             position = position_jitter(width = 0.25)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  theme_classic() +
  labs(x = "Gene",
       y = "Δ allele frequency",
       size = "CADD",
       colour = "Impact",
       shape = "Direction")

# Save

ggsave("delta_cadd_plot.png", p_cadd, width = 10, height = 8, dpi = 300)

# =========================
# END
# =========================
