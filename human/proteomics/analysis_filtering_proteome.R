# Title: Image preparation for proteomics data based on GS filtering
# Author: Jayme Cohen
#16/10/2025

library(readr)
library(dplyr)
library(pheatmap)
library(gridExtra)
library(VennDiagram)
library(ggseqlogo)
library(ggrepel)
library(data.table)
library(tidyverse)
library(biomaRt)
library(readr)
library(openxlsx)
library(knitr)
library(kableExtra)
library(tidyr)
library(ggplot2)
library(patchwork)
library(stringr)

GS_dat <- readxl::read_excel("GS_dat.xlsx", skip = 2)
EXP1 <- readxl::read_excel("EXP1_2025.xlsx")
EXP2 <- readxl::read_excel("EXP2_2025.xlsx")
EXP3 <- readxl::read_excel("EXP3_2025_correct.xlsx")

#rename the columns in the GS table for clarity
GS_dat <- GS_dat %>%
  rename(
    Donor1_ratio = `B/A...3`,
    Donor2_ratio = `D/C...4`,
    Donor3_ratio = `2vs1...6`,
    Donor4_ratio = `6vs5...7`,
    Donor5_ratio = `Exp3=2/1`,
    P_val_don1 = `B/A...9`,
    P_val_don2 = `D/C...10`,
    P_val_don3 = `2vs1...12`,
    P_val_don4 = `6vs5...13`,
    P_val_don5 = 'Exp3=B/A', 
    EXP1_Unique_pep = Exp1,
    EXP2_Unique_pep = Exp2,
    EXP3_Unique_pep = Exp3
  )

# Helper function to simplify colname cleanup
clean_ab_colnames <- function(df, mapping) {
  names(df) <- str_replace_all(names(df), mapping)
  return(df)
}

# EXP1
exp1_cols <- EXP1 %>%
  dplyr::select(
    Accession,
    `Abundance: F16: 127N, Sample, A_don1pre`,
    `Abundance: F16: 127C, Sample, B_don1post`,
    `Abundance: F16: 128N, Sample, C_don2pre`,
    `Abundance: F16: 128C, Sample, D_don2post`,
    `Abundances (Normalized): F16: 127N, Sample, A_don1pre`,
    `Abundances (Normalized): F16: 127C, Sample, B_don1post`,
    `Abundances (Normalized): F16: 128N, Sample, C_don2pre`,
    `Abundances (Normalized): F16: 128C, Sample, D_don2post`
  ) %>%
  rename(
    Donor1_pre = `Abundance: F16: 127N, Sample, A_don1pre`,
    Donor1_post = `Abundance: F16: 127C, Sample, B_don1post`,
    Donor1_norm_pre = `Abundances (Normalized): F16: 127N, Sample, A_don1pre`,
    Donor1_norm_post = `Abundances (Normalized): F16: 127C, Sample, B_don1post`,
    Donor2_pre = `Abundance: F16: 128N, Sample, C_don2pre`,
    Donor2_post = `Abundance: F16: 128C, Sample, D_don2post`,
    Donor2_norm_pre = `Abundances (Normalized): F16: 128N, Sample, C_don2pre`,
    Donor2_norm_post = `Abundances (Normalized): F16: 128C, Sample, D_don2post`
  )

# EXP2
exp2_cols <- EXP2 %>%
  dplyr::select(
    Accession,
    "Abundance: F15: 127N, Sample, D1pre",
    "Abundance: F15: 127C, Sample, D1post", 
    "Abundance: F15: 129N, Sample, D3pre",
    `Abundance: F15: 129C, Sample, D3post`,
    "Abundances (Normalized): F15: 127N, Sample, D1pre",
    "Abundances (Normalized): F15: 127C, Sample, D1post", 
   "Abundances (Normalized): F15: 129N, Sample, D3pre",
    `Abundances (Normalized): F15: 129C, Sample, D3post`
  ) %>%
  rename(
    Donor3_pre = "Abundance: F15: 127N, Sample, D1pre",
    Donor3_post = "Abundance: F15: 127C, Sample, D1post", 
    Donor4_pre = "Abundance: F15: 129N, Sample, D3pre",
   Donor4_post = `Abundance: F15: 129C, Sample, D3post`,
    Donor3_norm_pre = "Abundances (Normalized): F15: 127N, Sample, D1pre",
    Donor3_norm_post = "Abundances (Normalized): F15: 127C, Sample, D1post", 
    Donor4_norm_pre = "Abundances (Normalized): F15: 129N, Sample, D3pre",
    Donor4_norm_post = `Abundances (Normalized): F15: 129C, Sample, D3post`
  )

# EXP3
exp3_cols <- EXP3 %>%
  dplyr::select(
    Accession,
    "Abundance: F3: 127N, Sample, 1A",
    "Abundance: F3: 127C, Sample, 1B",
    "Abundances (Normalized): F3: 127N, Sample, 1A",
    "Abundances (Normalized): F3: 127C, Sample, 1B"
  ) %>%
  rename(
    Donor5_pre = "Abundance: F3: 127N, Sample, 1A",
    Donor5_post = "Abundance: F3: 127C, Sample, 1B",
    Donor5_norm_pre = "Abundances (Normalized): F3: 127N, Sample, 1A",
    Donor5_norm_post = "Abundances (Normalized): F3: 127C, Sample, 1B"
  )

Full_data <- GS_dat %>%
  left_join(exp1_cols, by = "Accession") %>%
  left_join(exp2_cols, by = "Accession") %>%
  left_join(exp3_cols, by = "Accession")

write.csv(Full_data, "Full_Proteomics_Data_unfiltered.csv", row.names = FALSE)

Filtered_data <- Full_data %>%
  filter(
    (`EXP1_Unique_pep` >= 2 | is.na(`EXP1_Unique_pep`)) &
      (`EXP2_Unique_pep` >= 2 | is.na(`EXP2_Unique_pep`)) &
      (`EXP3_Unique_pep` >= 2 | is.na(`EXP3_Unique_pep`)) &
      `#present in` >= 3
  ) #filter to remove below two unique peptides, keeping proteins that dont show up in all exp

Filtered_data <- Filtered_data %>%
  distinct(Accession, .keep_all = TRUE) #remove any duplicates 

GS_data_filtered <- GS_dat %>%
  filter(
    (`EXP1_Unique_pep` >= 2 | is.na(`EXP1_Unique_pep`)) &
      (`EXP2_Unique_pep` >= 2 | is.na(`EXP2_Unique_pep`)) &
      (`EXP3_Unique_pep` >= 2 | is.na(`EXP3_Unique_pep`)) &
      `#present in` >= 3
  ) #check against the original filter performed by GS

setdiff(Filtered_data$Accession, GS_data_filtered$Accession)
setdiff(GS_data_filtered$Accession, Filtered_data$Accession)
#should both be zero- YES!

write.csv(Filtered_data, "Full_Proteomics_Data_filtered.csv", row.names = FALSE)

save.image(file = "rdata_filtered.RData")
#___________________________________________________________________________

#Experiment subsets and single donor subsets
metadata_cols <- c("Accession", "Description", "Median", "#present in", 
                   "#Sign_P", "#down", "#up", "#empty")

donor1_2_cols <- c(
  "Donor1_ratio", "P_val_don1", "Donor1_pre", "Donor1_post",
  "Donor1_norm_pre", "Donor1_norm_post",
  "Donor2_ratio", "P_val_don2", "Donor2_pre", "Donor2_post",
  "Donor2_norm_pre", "Donor2_norm_post"
)

EXP1_final <- Filtered_data %>%
  dplyr::select(all_of(c(metadata_cols, donor1_2_cols)))

donor3_4_cols <- c(
  "Donor3_ratio", "P_val_don3", "Donor3_pre", "Donor3_post",
  "Donor3_norm_pre", "Donor3_norm_post",
  "Donor4_ratio", "P_val_don4", "Donor4_pre", "Donor4_post",
  "Donor4_norm_pre", "Donor4_norm_post"
)

EXP2_final <- Filtered_data %>%
  dplyr::select(all_of(c(metadata_cols, donor3_4_cols)))

donor5_cols <- c(
  "Donor5_ratio", "P_val_don5", "Donor5_pre", "Donor5_post",
  "Donor5_norm_pre", "Donor5_norm_post"
)

EXP3_final <- Filtered_data %>%
  dplyr::select(all_of(c(metadata_cols, donor5_cols)))

donor1_cols <- c("Donor1_ratio", "P_val_don1", "Donor1_pre", "Donor1_post",
                 "Donor1_norm_pre", "Donor1_norm_post")
Donor1_Final <- Filtered_data %>%
  dplyr::select(all_of(c(metadata_cols, donor1_cols)))

donor2_cols <- c("Donor2_ratio", "P_val_don2", "Donor2_pre", "Donor2_post",
                 "Donor2_norm_pre", "Donor2_norm_post")

Donor2_Final <- Filtered_data %>%
  dplyr::select(all_of(c(metadata_cols, donor2_cols)))

Donor3_cols <- c("Donor3_ratio", "P_val_don3", "Donor3_pre", "Donor3_post",
                 "Donor3_norm_pre", "Donor3_norm_post")
Donor3_Final <- Filtered_data %>%
  dplyr::select(all_of(c(metadata_cols, Donor3_cols)))

donor4_cols <- c("Donor4_ratio", "P_val_don4", "Donor4_pre", "Donor4_post",
                 "Donor4_norm_pre", "Donor4_norm_post")

Donor4_Final <- Filtered_data %>%
  dplyr::select(all_of(c(metadata_cols, donor4_cols)))

Donor1_Final <- Donor1_Final %>%
  filter(!is.na(Donor1_ratio))

Donor2_Final <- Donor2_Final %>%
  filter(!is.na(Donor2_ratio))

Donor3_Final <- Donor3_Final %>%
  filter(!is.na(Donor3_ratio))

Donor4_Final <- Donor4_Final %>%
  filter(!is.na(Donor4_ratio))

Donor5_Final <- EXP3_final %>%
  filter(!is.na(Donor5_ratio))
#___________________________________________________________________________
#make overall lists of proteins up and down in the selected for (post) group
#significant in at least 1/5, same trend in 3/5 minimum.

overall_up <- Filtered_data %>%
  filter(`#Sign_P` >= 1, `#up` >= 3)

overall_down <- Filtered_data %>%
  filter(`#Sign_P` >= 1, `#down` >= 3)

#Shortlist from the overall- cutoff of median ratio of 0.2 and 2 shared by 3/5

Filtered_data <- Filtered_data %>%
  mutate(
    num_strong_up = rowSums(across(contains("_ratio"), ~ . >= 2, .names = "up"), na.rm = TRUE),
    num_strong_down = rowSums(across(contains("_ratio"), ~ . <= 0.2, .names = "down"), na.rm = TRUE)
  )

Shortlist_up <- Filtered_data %>%
  filter(`#Sign_P` >= 1, num_strong_up >= 3)

Shortlist_down <- Filtered_data %>%
  filter(`#Sign_P` >= 1, num_strong_down >= 3)

Shortlist_down <- Shortlist_down %>%
  mutate(Direction = "Down")

Shortlist_up <- Shortlist_up %>%
  mutate(Direction = "Up")

Shortlist_overall <- bind_rows(Shortlist_down, Shortlist_up)

write.csv(Shortlist_overall, "Proteomics_shortlist.csv", row.names = FALSE)
#__________________________________________________________________________
#Up and down per experiment, significant in at least 1/2 and ratio >1 or <1 in both

exp1_up <- EXP1_final %>%
  filter(
    Donor1_ratio > 1, Donor2_ratio > 1,
    P_val_don1 < 0.05 | P_val_don2 < 0.05
  )

exp1_down <- EXP1_final %>%
  filter(
    Donor1_ratio < 1, Donor2_ratio < 1,
    P_val_don1 < 0.05 | P_val_don2 < 0.05
  )

exp2_up <- EXP2_final %>%
  filter(
    Donor3_ratio > 1, Donor4_ratio > 1,
    P_val_don3 < 0.05 | P_val_don4 < 0.05
  )

exp2_down <- EXP2_final %>%
  filter(
    Donor3_ratio < 1, Donor4_ratio < 1,
    P_val_don3 < 0.05 | P_val_don4 < 0.05
  )

exp3_up <- EXP3_final %>%
  filter(Donor5_ratio > 1, P_val_don5 < 0.05)

exp3_down <- EXP3_final %>%
  filter(Donor5_ratio < 1, P_val_don5 < 0.05)
#__________________________________________________________________________
#per donor must be significant and meet ratio

donor1_up <- Donor1_Final %>%
  filter(Donor1_ratio > 1, P_val_don1 < 0.05)

donor1_down <- Donor1_Final %>%
  filter(Donor1_ratio < 1, P_val_don1 < 0.05)

donor2_up <- Donor2_Final %>%
  filter(Donor2_ratio > 1, P_val_don2 < 0.05)

donor2_down <- Donor2_Final %>%
  filter(Donor2_ratio < 1, P_val_don2 < 0.05)

donor3_up <- Donor3_Final %>%
  filter(Donor3_ratio > 1, P_val_don3 < 0.05)

donor3_down <- Donor3_Final %>%
  filter(Donor3_ratio < 1, P_val_don3 < 0.05)

donor4_up <- Donor4_Final %>%
  filter(Donor4_ratio > 1, P_val_don4 < 0.05)

donor4_down <- Donor4_Final %>%
  filter(Donor4_ratio < 1, P_val_don4 < 0.05)

donor5_up <- Donor5_Final %>%
  filter(Donor5_ratio > 1, P_val_don5 < 0.05)

donor5_down <- Donor5_Final %>%
  filter(Donor5_ratio < 1, P_val_don5 < 0.05)


#___________________________________________________________________________

#make top5 based on whole dataset- significant in at least 4/5,
#present in all, same trend in all, Median less than or equal to 0.2

top5_down <- Shortlist_down %>%
  filter(
    `#Sign_P` >= 4,
    `#present in` == 5,
    Median <= 0.2
  ) %>%
  arrange(Median) %>%         # sort smallest to largest median ratio
  slice_head(n = 5)

top5_up <- Shortlist_up %>%
  filter(
    `#Sign_P` >= 4,
    `#present in` == 5,
    Median >= 2
  ) %>%
  arrange(desc(Median)) %>%         # sort smallest to largest median ratio
  slice_head(n = 5)

# EXP1
EXP1_combine <- EXP1 %>%
  dplyr::mutate(across(
    dplyr::matches("^Abundance"), log2,  # only log2-transform numeric abundance columns
    .names = "{.col}"
  )) %>%
  dplyr::select(
    Accession, Description,
    dplyr::starts_with("Abundance"),
    dplyr::starts_with("P-Value"),
    dplyr::starts_with("Adj. P-Value")
  )

# EXP2
EXP2_combine <- EXP2 %>%
  dplyr::mutate(across(
    dplyr::matches("^Abundance"), log2,
    .names = "{.col}"
  )) %>%
  dplyr::select(
    Accession, Description,
    dplyr::starts_with("Abundance"),
    dplyr::starts_with("P-Value"),
    dplyr::starts_with("Adj. P-Value")
  )

# EXP3
EXP3_combine <- EXP3 %>%
  dplyr::mutate(across(
    dplyr::matches("^Abundance"), log2,
    .names = "{.col}"
  )) %>%
  dplyr::select(
    Accession, Description,
    dplyr::starts_with("Abundance"),
    dplyr::starts_with("P-Value"),
    dplyr::starts_with("Adj. P-Value")
  )


EXP1_fc <- EXP1_combine %>%
  group_by(Accession) %>%
  mutate(log_fc1 = `Abundances (Normalized): F16: 127N, Sample, A_don1pre` - `Abundances (Normalized): F16: 127C, Sample, B_don1post`,
         log_fc2 = `Abundances (Normalized): F16: 128N, Sample, C_don2pre` - `Abundances (Normalized): F16: 128C, Sample, D_don2post`,
         log_pval1 = -1*log10(`P-Value: (B_don1post) / (A_don1pre)`),
         log_pval2 = -1*log10(`P-Value: (D_don2post) / (C_don2pre)`))

EXP2_fc <- EXP2_combine %>%
  group_by(Accession) %>%
  mutate(log_fc1 = `Abundances (Normalized): F15: 127N, Sample, D1pre` - `Abundances (Normalized): F15: 127C, Sample, D1post`,
         log_fc2 = `Abundances (Normalized): F15: 129N, Sample, D3pre` - `Abundances (Normalized): F15: 129C, Sample, D3post`,
         log_pval1 = -1*log10(`P-Value: (D1post) / (D1pre)`),
         log_pval2 = -1*log10(`P-Value: (D3post) / (D3pre)`))

EXP3_fc <- EXP3 %>%
  dplyr::group_by(Accession) %>%
  dplyr::mutate(
    log_fc = log2(
      (
        (`Abundances (Normalized): F3: 127N, Sample, 1A` +
           `Abundances (Normalized): F4: 127N, Sample, 1A`) / 2
      ) /
        (
          (`Abundances (Normalized): F3: 127C, Sample, 1B` +
             `Abundances (Normalized): F4: 127C, Sample, 1B`) / 2
        )
    ),
    log_pval = -log10(`Abundance Ratio P-Value: (1B) / (1A)`)
  )





# Top 5 strongly less abundant after selection (smallest Median_Ratio)
top5_down <- shortlist_down %>%
  arrange(desc(num_significant), Median_Ratio) %>%  # Sort by increasing Median_Ratio
  slice_head(n = 5)

top5_up <- shortlist_up %>%
  arrange(num_significant, Median_Ratio) %>%  # Sort by increasing Median_Ratio
  slice_head(n = 5)

save.image(file = "rdata_filtered_shortlists_finished.RData")

#==============================#
# Step 6: plot the volcanos
#==============================#
# Make sure merged_filtered has short_name
merged_filtered <- Filtered_data %>%
  mutate(
    short_name = ifelse(
      !is.na(Description) & str_detect(Description, " OS="),
      str_remove(str_extract(Description, "^.*? OS="), " OS=$"),
      Description
    )
  )


# Attach labels to top5 lists
top5_down <- Shortlist_down %>%
  filter(`#present in` == 5) %>%   # ✅ require presence in all 5 donors
  arrange(desc(`#Sign_P`), `#present in`, Median) %>%
  slice_head(n = 5)

top5_down <- top5_down %>%
  left_join(merged_filtered %>% dplyr::select(Accession, short_name), by = "Accession") %>%
  mutate(label = paste0(Accession, " (", short_name, ")"))


top5_up <- Shortlist_up %>%
  arrange(`#Sign_P`, Median) %>%  # Sort by increasing Median_Ratio
  slice_head(n = 5)

top5_up <- top5_up %>%
  left_join(merged_filtered %>% dplyr::select(Accession, short_name), by = "Accession") %>%
  mutate(label = paste0(Accession, " (", short_name, ")"))

shortlist <- bind_rows(top5_down, top5_up) %>%
  distinct(Accession, .keep_all = TRUE)

highlight_accessions <- shortlist$Accession


# Generate donor-specific data frames with flipped FC
df1 <- Donor1_Final %>%
  mutate(
    log_fc1 = -log2(Donor1_ratio),
    log_pval1 = -log10(P_val_don1)
  )

df2 <- Donor2_Final %>%
  mutate(
    log_fc2 = -log2(Donor2_ratio),
    log_pval2 = -log10(P_val_don2)
  )

df3 <- Donor3_Final %>%
  mutate(
    log_fc1 = -log2(Donor3_ratio),
    log_pval1 = -log10(P_val_don3)
  )

df4 <- Donor4_Final %>%
  mutate(
    log_fc2 = -log2(Donor4_ratio),
    log_pval2 = -log10(P_val_don4)
  )

df5 <- Donor5_Final %>%
  mutate(
    log_fc = -log2(Donor5_ratio),
    log_pval = -log10(P_val_don5)
  )

# Function to prepare a single volcano dataset
prepare_volcano_data <- function(df, logfc_col, pval_col) {
  df <- df %>%
    left_join(merged_filtered %>% dplyr::select(Accession, short_name), by = "Accession") %>%
    mutate(
      log2FoldChange = -!!sym(logfc_col),
      pval = 10^(-!!sym(pval_col)),                # convert -log10 back to p for consistency
      neg_log10_pval = !!sym(pval_col),
      significance = case_when(
        neg_log10_pval > -log10(0.05) & log2FoldChange >= 0.5  ~ "Higher Abundance in T4",
        neg_log10_pval > -log10(0.05) & log2FoldChange <= -0.5 ~ "Lower Abundance in T4",
        TRUE ~ "Not Significant"
      ),
      label = ifelse(
        Accession %in% highlight_accessions,
        paste0(Accession, " (", short_name, ")"),
        NA_character_
      )
    )
  
  return(df)
}

# Unified plotting function with proper formatting
plot_volcano_unified <- function(df, title, show_legend = FALSE) {
  ggplot(df, aes(x = log2FoldChange, y = neg_log10_pval)) +
    geom_point(
      aes(color = significance),
      alpha = 0.6, size = 2.5,           # slightly larger points
      position = position_jitter(width = 0.15, height = 0.15)
    ) +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
    scale_color_manual(
      values = c(
        "Higher Abundance in T4" = "purple",
        "Lower Abundance in T4" = "darkorange1",
        "Not Significant" = "steelblue2"
      ),
      guide = if (show_legend) "legend" else "none"
    ) +
    geom_text_repel(
      data = filter(df, !is.na(label)),
      aes(label = label),
      size = 4,                # 👈 THIS controls the label text size
      fontface = "bold",
      box.padding = 0.35,
      point.padding = 0.25,
      force = 6,
      max.iter = 5000,
      max.overlaps = Inf,
      segment.size = 0.25,
      segment.alpha = 0.4,
      min.segment.length = 0
    ) +
    labs(
      title = title,
      x = "-log2 Fold Change",
      y = "-log10 p-value"
    ) +
    coord_cartesian(clip = "off") +
    theme_minimal(base_size = 50) +      # <-- this globally increases font sizes
    theme(
      legend.position = if (show_legend) "right" else "none",
      legend.text = element_text(size = 13),
      legend.title = element_text(size = 14, face = "bold"),
      axis.title = element_text(size = 15),
      axis.text = element_text(size = 13),
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      panel.grid.minor = element_blank(),
      axis.title.x = element_text(face = "bold"),
      axis.title.y = element_text(face = "bold")
    )
}

# Prepare the five data frames
plot_data1 <- prepare_volcano_data(df1, "log_fc1", "log_pval1")
plot_data2 <- prepare_volcano_data(df2, "log_fc2", "log_pval2")
plot_data3 <- prepare_volcano_data(df3, "log_fc1", "log_pval1")
plot_data4 <- prepare_volcano_data(df4, "log_fc2", "log_pval2")
plot_data5 <- prepare_volcano_data(df5, "log_fc",  "log_pval")

# Create the five plots
p1 <- plot_volcano_unified(plot_data1, "Methyl Cellulose Donor 1")
p2 <- plot_volcano_unified(plot_data2, "Methyl Cellulose Donor 2")
p3 <- plot_volcano_unified(plot_data3, "Swim-up Donor 3")
p4 <- plot_volcano_unified(plot_data4, "Swim-up Donor 4")
p5 <- plot_volcano_unified(plot_data5, "Swim-up Run 2 Donor 3")


print(p1)
print(p2)
print(p3)
print(p4)
print(p5)

p1 <- p1 + ggtitle("Donor 1")
p2 <- p2 + ggtitle("Donor 2")
p3 <- p3 + ggtitle("Donor 3")
p4 <- p4 + ggtitle("Donor 4")
p5 <- p5 + ggtitle("Run 2 Donor 3")

plots <- list(p1, p2, p3, p4, p5)
names(plots) <- c("p1_MethylCellulose_Donor1",
                  "p2_MethylCellulose_Donor2",
                  "p3_SwimUp_Donor3",
                  "p4_SwimUp_Donor4",
                  "p5_SwimUpRepeat_Donor3")

for (n in names(plots)) {
  plots[[n]] <- plots[[n]] + theme(
    text = element_text(size = 20),
    axis.title = element_text(size = 30, face = "bold"),
    axis.text  = element_text(size = 30),
    plot.title = element_text(size = 30, face = "bold", hjust = 0.5)
  )
  
  ggsave(
    filename = paste0(n, "v.tiff"),
    plot = plots[[n]],
    width = 13,
    height = 12,
    units = "in",
    scale = 0.8,
    dpi = 600,
    bg = "white",
    compression = "lzw"
  )
}


#___________________________________________________________________________
#Add metadata for GO terms and location annotations 

library(biomaRt)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(dplyr)

all_accessions <- unique(Filtered_data$Accession)

mart <-mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

mapping <- getBM(
  attributes = c("uniprotswissprot", "ensembl_gene_id", "entrezgene_id"),
  filters = "uniprotswissprot",
  values = all_accessions,
  mart = mart
)

mapping <- mapping %>%
  filter(!is.na(entrezgene_id)) %>%
  distinct(uniprotswissprot, .keep_all = TRUE)

gene_entrez <- unique(mapping$entrezgene_id)

library(UniProt.ws)

up <- UniProt.ws(taxId = 9606)  # 9606 = human
columns(up)
keytypes(up)

annotations <- select(
  up,
  keys = all_accessions,
  columns = c("cc_subcellular_location"),
  keytype = "UniProtKB"
)

library(dplyr)
library(stringr)

cleaned_annotations <- annotations %>%
  mutate(
    Clean_Location = str_remove(Subcellular.location..CC., "^SUBCELLULAR LOCATION: "),
    Clean_Location = str_remove(Clean_Location, "\\{.*$"),
    Clean_Location = str_trim(Clean_Location)
  )

Shortlist_overall_annotated <- Shortlist_overall %>%
  left_join(
    cleaned_annotations %>%
      dplyr::select(From, Clean_Location),
    by = c("Accession" = "From")
  )

ranked_shortlist <- Shortlist_overall_annotated %>%
  mutate(
    Membrane_flag = if_else(str_detect(tolower(Clean_Location), "membrane"), 1, 0)
  ) %>%
  arrange(desc(Membrane_flag), Median)

write.csv(ranked_shortlist, "Shortlist_ranked_annotated.csv", row.names = FALSE)

Filtered_data_annotated <- Filtered_data %>%
  left_join(
    cleaned_annotations %>%
      dplyr::select(From, Clean_Location),
    by = c("Accession" = "From")
  )
write.csv(Filtered_data_annotated, "Filtered_proteome_annotated.csv", row.names = FALSE)

#___________________________________________________________________________
#Make GO term plots per EXP up and EXP down 
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(patchwork)
library(dplyr)

run_enrichment_for_plotting <- function(accession_list, mapping_df) {
  gene_entrez <- mapping_df %>%
    filter(uniprotswissprot %in% accession_list) %>%
    pull(entrezgene_id) %>%
    unique() %>%
    na.omit()
  
  ego_bp <- enrichGO(gene = gene_entrez, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
                     ont = "BP", pvalueCutoff = 0.05, readable = TRUE)
  
  ego_cc <- enrichGO(gene = gene_entrez, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
                     ont = "CC", pvalueCutoff = 0.05, readable = TRUE)
  
  ego_mf <- enrichGO(gene = gene_entrez, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
                     ont = "MF", pvalueCutoff = 0.05, readable = TRUE)
  
  ekegg <- enrichKEGG(gene = gene_entrez, organism = "hsa", pvalueCutoff = 0.05)
  
  list(
    BP = ego_bp,
    CC = ego_cc,
    MF = ego_mf,
    KEGG = ekegg
  )
}

res_exp1_up <- run_enrichment_for_plotting(exp1_up$Accession, mapping)
res_exp2_up <- run_enrichment_for_plotting(exp2_up$Accession, mapping)
res_exp3_up <- run_enrichment_for_plotting(exp3_up$Accession, mapping)

library(patchwork)
library(enrichplot)

# Combine GO:BP
bp_plot <- (
  dotplot(res_exp1_up$BP, showCategory = 10, title = "Methyl Cellulose") |
    dotplot(res_exp2_up$BP, showCategory = 10, title = "Swim Up") |
    dotplot(res_exp3_up$BP, showCategory = 10, title = "Swim Up TMT #2")
)+
  plot_annotation(title = "GO: Biological Processes Enrichment – Proteins More Abundant in Selected Sperm") &
  theme(
    text = element_text(size = 30),          # ⬅️ increases *all* text
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 14)
  )

print(bp_plot)
ggsave(
  "GO_BP_up_largefonts.png",
  plot = bp_plot,
  width = 8, height = 4, units = "in",
  dpi = 600,
  scale = 2.5 # ⬅️ ensures fonts stay the same regardless of image size
)
# Combine GO:CC
cc_plot <- (
  dotplot(res_exp1_up$CC, showCategory = 10, title = "Methyl Cellulose") |
    dotplot(res_exp2_up$CC, showCategory = 10, title = "Swim Up") |
    dotplot(res_exp3_up$CC, showCategory = 10, title = "Swim Up TMT #2")
)+
  plot_annotation(title = "GO: Cellular Component Enrichment – Proteins More Abundant in Selected Sperm")&
  theme(
    text = element_text(size = 30),          # ⬅️ increases *all* text
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 14)
  )
print(cc_plot)
ggsave(
  "GO_CC_up_largefonts.png",
  plot = cc_plot,
  width = 8, height = 4, units = "in",
  dpi = 600,
  scale = 2.5        # ⬅️ ensures fonts stay the same regardless of image size
)

# Combine GO:MF
mf_plot <- (
  dotplot(res_exp1_up$MF, showCategory = 10, title = "Methyl Cellulose") |
    dotplot(res_exp2_up$MF, showCategory = 10, title = "Swim Up") |
    dotplot(res_exp3_up$MF, showCategory = 10, title = "Swim Up TMT #2")
)+
  plot_annotation(title = "GO: Molecular Function Enrichment – Proteins More Abundant in Selected Sperm")&
  theme(
    text = element_text(size = 30),          # ⬅️ increases *all* text
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 14)
  )
print(mf_plot)
ggsave(
  "GO_MF_up_largefonts.png",
  plot = mf_plot,
  width = 8, height = 4.5, units = "in",
  dpi = 600,
  scale = 2.6        # ⬅️ ensures fonts stay the same regardless of image size
)


# Combine KEGG
kegg_plot <- (
  dotplot(res_exp1_up$KEGG, showCategory = 10, title = "Methyl Cellulose") |
    dotplot(res_exp2_up$KEGG, showCategory = 10, title = "Swim Up") |
    dotplot(res_exp3_up$KEGG, showCategory = 10, title = "Swim Up TMT #2")
)+
  plot_annotation(title = "KEGG Pathway Enrichment – Proteins More Abundant in Selected Sperm")&
  theme(
    text = element_text(size = 30),          # ⬅️ increases *all* text
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 14)
  )
print(kegg_plot)
ggsave(
  "GO_KEGG_up_largefonts.png",
  plot = kegg_plot,
  width = 8, height = 4, units = "in",
  dpi = 600,
  scale = 2.65       # ⬅️ ensures fonts stay the same regardless of image size
)



# Run enrichment for downregulated accessions
res_exp1_down <- run_enrichment_for_plotting(exp1_down$Accession, mapping)
res_exp2_down <- run_enrichment_for_plotting(exp2_down$Accession, mapping)
res_exp3_down <- run_enrichment_for_plotting(exp3_down$Accession, mapping)

bp_plot_down <- (
  dotplot(res_exp1_down$BP, showCategory = 10, title = "Methyl Cellulose") |
    dotplot(res_exp2_down$BP, showCategory = 10, title = "Swim Up") |
    dotplot(res_exp3_down$BP, showCategory = 10, title = "Swim Up TMT #2")
) +
  plot_annotation(title = "GO: Biological Processes Enrichment – Proteins Less Abundant in Selected Sperm")&
  theme(
    text = element_text(size = 30),          # ⬅️ increases *all* text
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 14)
  )
print(bp_plot_down)
ggsave(
  "GO_BP_down_largefonts.png",
  plot = bp_plot_down,
  width = 14, height = 6, units = "in",
  dpi = 600,
  scale = 1.35        # ⬅️ ensures fonts stay the same regardless of image size
)

cc_plot_down <- (
  dotplot(res_exp1_down$CC, showCategory = 10, title = "Methyl Cellulose") |
    dotplot(res_exp2_down$CC, showCategory = 10, title = "Swim Up") |
    dotplot(res_exp3_down$CC, showCategory = 10, title = "Swim Up TMT #2")
) +
  plot_annotation(title = "GO: Cellular Component Enrichment – Proteins Less Abundant in Selected Sperm")&
  theme(
    text = element_text(size = 30),          # ⬅️ increases *all* text
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 14)
  )
print(CC_plot_down)
ggsave(
  "GO_CC_down_largefonts.png",
  plot = cc_plot_down,
  width = 14, height = 6, units = "in",
  dpi = 600,
  scale = 1.45        # ⬅️ ensures fonts stay the same regardless of image size
)

mf_plot_down <- (
  dotplot(res_exp1_down$MF, showCategory = 10, title = "Methyl Cellulose") |
    dotplot(res_exp2_down$MF, showCategory = 10, title = "Swim Up") |
    dotplot(res_exp3_down$MF, showCategory = 10, title = "Swim Up TMT #2")
) +
  plot_annotation(title = "GO: Molecular Function Enrichment – Proteins Less Abundant in Selected Sperm")&
  theme(
    text = element_text(size = 30),          # ⬅️ increases *all* text
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 14)
  )
print(mf_plot_down)
ggsave(
  "GO_MF_down_largefonts.png",
  plot = mf_plot_down,
  width = 14, height = 6, units = "in",
  dpi = 600,
  scale = 1.45        # ⬅️ ensures fonts stay the same regardless of image size
)

kegg_plot_down <- (
  dotplot(res_exp1_down$KEGG, showCategory = 10, title = "Methyl Cellulose") |
    dotplot(res_exp2_down$KEGG, showCategory = 10, title = "Swim Up") |
    dotplot(res_exp3_down$KEGG, showCategory = 10, title = "Swim Up TMT #2")
) +
  plot_annotation(title = "KEGG Pathway Enrichment – Proteins Less Abundant in Selected Sperm")&
  theme(
    text = element_text(size = 30),          # ⬅️ increases *all* text
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 14)
  )
print(kegg_plot_down)
ggsave(
  "GO_kegg_down_largefonts.png",
  plot = kegg_plot_down,
  width = 14, height = 6, units = "in",
  dpi = 600,
  scale = 1.45        # ⬅️ ensures fonts stay the same regardless of image size
)

res_shortlist_down <- run_enrichment_for_plotting(Shortlist_down$Accession, mapping)
res_shortlist_up <- run_enrichment_for_plotting(Shortlist_up$Accession, mapping)

shortlist_up_plot <- (
  dotplot(res_shortlist_up$BP, showCategory = 10, title = "GO: Biological Process") |
    dotplot(res_shortlist_up$CC, showCategory = 10, title = "GO: Cellular Component")
) /
  (
    dotplot(res_shortlist_up$MF, showCategory = 10, title = "GO: Molecular Function") |
      dotplot(res_shortlist_up$KEGG, showCategory = 10, title = "KEGG Pathways")
  ) +
  plot_annotation(title = "Functional Enrichment – Proteins More Abundant in Selected Sperm (Shortlist)")

library(ggplot2)

# Safe fallback dotplot wrapper
safe_dotplot <- function(enrichment_obj, title) {
  if (is.null(enrichment_obj) || nrow(as.data.frame(enrichment_obj)) == 0) {
    ggplot() + 
      annotate("text", x = 0.5, y = 0.5, label = "No enriched terms", size = 5, hjust = 0.5) +
      ggtitle(title) +
      theme_void()
  } else {
    dotplot(enrichment_obj, showCategory = 10, title = title)
  }
}

# Shortlist UP combined plot
shortlist_up_plot <- (
  safe_dotplot(res_shortlist_up$BP, "GO: Biological Process") |
    safe_dotplot(res_shortlist_up$CC, "GO: Cellular Component")
) /
  (
    safe_dotplot(res_shortlist_up$MF, "GO: Molecular Function") |
      safe_dotplot(res_shortlist_up$KEGG, "KEGG Pathways")
  ) +
  plot_annotation(title = "Functional Enrichment – Proteins More Abundant in Selected Sperm (Shortlist)")&
  theme(
    text = element_text(size = 30)         # ⬅️ increases *all* text
  )
print(shortlist_up_plot)
ggsave(
  "GO_shortlist_up_largefonts.png",
  plot = shortlist_up_plot,
  width = 14, height = 7, units = "in",
  dpi = 600,
  scale = 1.4        # ⬅️ ensures fonts stay the same regardless of image size
)
#down
shortlist_down_plot <- (
  safe_dotplot(res_shortlist_down$BP, "GO: Biological Process") |
    safe_dotplot(res_shortlist_down$CC, "GO: Cellular Component")
) /
  (
    safe_dotplot(res_shortlist_down$MF, "GO: Molecular Function") |
      safe_dotplot(res_shortlist_down$KEGG, "KEGG Pathways")
  ) +
  plot_annotation(title = "Functional Enrichment – Proteins Less Abundant in Selected Sperm (Shortlist)")&
  theme(
    text = element_text(size = 30)         # ⬅️ increases *all* text
  )
print(shortlist_down_plot)
ggsave(
  "GO_shortlist_down_largefonts.png",
  plot = shortlist_down_plot,
  width = 14, height = 8, units = "in",
  dpi = 600,
  scale = 1.4        # ⬅️ ensures fonts stay the same regardless of image size
)

#___________________________________________________________________________
#make Venn diagrams for just proteomics 

# Load required libraries
library(tidyverse)
library(VennDiagram)
library(gridExtra)
library(readxl)
library(biomaRt)

venn_donors <- list(
  "D1- Methyl Cellulose" = unique(Donor1_Final$Accession),
  "D2- Methyl Cellulose" = unique(Donor2_Final$Accession),
  "D3- Swim Up" = unique(Donor3_Final$Accession),
  "D4- Swim Up" = unique(Donor4_Final$Accession),
  "D3- Swim Up TMT#2" = unique(Donor5_Final$Accession)
)
venn_plot_all <- venn.diagram(
  venn_donors, NULL,
  col = "transparent",
  fill = c("cornflowerblue", "green", "yellow", "red", "violet"),
  alpha = 0.50,
  cex = 1.4,
  fontfamily = "sans",
  fontface = "bold",
  cat.col = c("darkblue", "darkgreen", "darkorange", "maroon", "purple"),
  cat.cex = 1.2,
  cat.fontfamily = "sans",
  main = "Proteins Identified per Donor",
  main.fontfamily = "sans",
  print.mode = "raw"
)
grid.arrange(gTree(children = venn_plot_all))


library(readr)
DM <- read_csv("Daniels Data GeneID .csv")

Filtered_data <- Filtered_data %>%
  left_join(mapping, by = c("Accession" = "uniprotswissprot"))

proteomics_ensembl <- Filtered_data %>%
  filter(!is.na(ensembl_gene_id)) %>%
  pull(ensembl_gene_id) %>%
  unique()

daniels_all <- unique(DM$Gene)

venn_proteo_vs_geno <- list(
  "Proteomics (Filtered)" = proteomics_ensembl,
  "Genomics" = daniels_all
)

DM_SW_only <- DM %>%
  filter(Experiment == "SW", !is.na(Gene)) %>%
  pull(Gene) %>%
  unique()

DM_MC_only <- DM %>%
  filter(Experiment == "MC", !is.na(Gene)) %>%
  pull(Gene) %>%
  unique()

venn_exp_vs_daniel <- list(
  "Proteomics MC" = unique(EXP1_final$Accession),
  "Proteomics SW" = unique(EXP2_final$Accession),
  "Proteomics SW TMT#2" = unique(EXP3_final$Accession),
  "Genomics MC" = DM_MC_only,
  "Genomics SW" = DM_SW_only
)

venn_plot_genomevproteome <- venn.diagram(
  venn_exp_vs_daniel, NULL,
  col = "transparent",
  fill = c("cornflowerblue", "green", "yellow", "red", "violet"),
  alpha = 0.50,
  cex = 1.4,
  fontfamily = "sans",
  fontface = "bold",
  cat.col = c("darkblue", "darkgreen", "darkorange", "maroon", "purple"),
  cat.cex = 1.2,
  cat.fontfamily = "sans",
  main = "Proteins Identified per Donor",
  main.fontfamily = "sans",
  print.mode = "raw"
)
grid.arrange(gTree(children = venn_plot_genomevproteome))


venn_daniels_genomics <- list(
  "Methyl cellulose" = DM_MC_only,
  "Swim-up" = DM_SW_only
)

genomic_only_plot <- venn.diagram(
  list("Genomics MC" = DM_MC_only, "Genomics SW" = DM_SW_only), NULL,
  col = "transparent",
  fill = c("darkblue", "darkmagenta"),
  alpha = 0.5,
  cex = 1.4,
  cat.cex = 1.2,
  main = "Genomics Only Overlap (Ensembl IDs)",
  print.mode = "raw"
)
grid.arrange(gTree(children = genomic_only_plot))

summary_counts <- data.frame(
  Experiment = c("Exp1", "Exp2", "Exp3"),
  Donors = c("1–2", "3–4", "5"),
  Up = c(nrow(exp1_up), nrow(exp2_up), nrow(exp3_up)),
  Down = c(nrow(exp1_down), nrow(exp2_down), nrow(exp3_down))
)

summary_counts_donors <- data.frame(
  Donor = paste0("Donor", 1:5),
  Up = c(nrow(donor1_up), nrow(donor2_up), nrow(donor3_up), nrow(donor4_up), nrow(donor5_up)),
  Down = c(nrow(donor1_down), nrow(donor2_down), nrow(donor3_down), nrow(donor4_down), nrow(donor5_down))
)

summary_counts
summary_counts_donors
