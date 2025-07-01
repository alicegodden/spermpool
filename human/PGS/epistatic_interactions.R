# PGS Score Change Scatter Plots (Dumbbell Plots), Cross-Phenotype Correlation Analysis,
# and Single Top Negative Correlation Scatter Plot

# --- INSTALL & LOAD NECESSARY LIBRARIES ---
# Ensure patchwork, Hmisc, corrplot, and tibble are installed
if (!requireNamespace("patchwork", quietly = TRUE)) {
  install.packages("patchwork")
}
if (!requireNamespace("Hmisc", quietly = TRUE)) {
  install.packages("Hmisc")
}
if (!requireNamespace("corrplot", quietly = TRUE)) {
  install.packages("corrplot")
}
if (!requireNamespace("tibble", quietly = TRUE)) { # Added tibble
  install.packages("tibble")
}
library(dplyr)      # For data manipulation
library(readr)      # For reading data (read_tsv)
library(ggplot2)    # For plotting
library(stringr)    # For string manipulation
library(tools)      # For file path manipulation (basename)
library(tidyr)      # For pivot_longer, pivot_wider
library(purrr)      # For reduce (to combine data frames efficiently)
library(patchwork)  # For combining plots
library(Hmisc)      # For rcorr (correlation matrix with p-values)
library(corrplot)   # For visualizing correlation matrix
library(tibble)     # For column_to_rownames

# --- SET WORKING DIRECTORY ---
# <<< IMPORTANT: Adjust this path to where your input files and phenocodes file are located >>>
setwd("~/PycharmProjects/polygenic/PGS_scores_test/SU/T0")  

# --- CONFIGURATION ---
# >>> IMPORTANT: Adjust these to match your actual file structure and desired behavior <<<

# Input files (raw score .tsv.gz files, assumed to be paired C, O, C, O...)
input_files <- c(
  "GRCh37_D1T0_merged.dedup_bcftools__CHR_nc_to_chr.nochr.aaf_scores.tsv.gz",
  "GRCh37_D1T4_merged.dedup_RG_bcftools__CHR_nc_to_chr.nochr.aaf_scores.tsv.gz",
  "GRCh37_D2T0_merged.dedup_bcftools__CHR_nc_to_chr.nochr.aaf_scores.tsv.gz",
  "GRCh37_D2T4_merged.dedup_bcftools__CHR_nc_to_chr.nochr.aaf_scores.tsv.gz",
  "GRCh37_D4T0_merged.dedup_bcftools__CHR_nc_to_chr.nochr.aaf_scores.tsv.gz",
  "GRCh37_D4T4_merged.dedup_bcftools__CHR_nc_to_chr.nochr.aaf_scores.tsv.gz",
  "GRCh37_D6T0a_merged.dedup_bcftools__CHR_nc_to_chr.nochr.aaf_scores.tsv.gz",
  "GRCh37_D6T4a_merged.dedup_bcftools__CHR_nc_to_chr.nochr.aaf_scores.tsv.gz",
  "GRCh37_D6T0b_merged.dedup_bcftools__CHR_nc_to_chr.nochr.aaf_scores.tsv.gz",
  "GRCh37_D6T4b_merged.dedup_bcftools__CHR_nc_to_chr.nochr.aaf_scores.tsv.gz"  
)

# Path to your phenocodes file (should contain 'filename' and 'description' columns)
phenocodes_path <- "phenocodes"  

# Output composite plot filename for Dumbbell plots
output_composite_dumbbell_plot_filename <- "PGS_Change_Dumbbell_Composite_PerPairTopPhenotypes.png"  

# Output filename for Correlation Heatmap
output_correlation_heatmap_filename <- "PGS_Change_Correlation_Heatmap_TopPhenotypes.png"

# Output filename for the single most negatively correlated scatter plot
output_single_scatter_plot_filename <- "PGS_Change_TopNegativeCorrelation_Scatter.png"

# Trait selection parameters for Dumbbell Plots (per-pair)
PERCENTILE_THRESHOLD <- 99 # Filter by top N percentile of absolute differences (now per pair)
MAX_TRAITS_TO_PLOT_DUMBBELL <- 15    # Cap for readability on individual dumbbell plots

# Trait selection parameters for Correlation Analysis (global)
MAX_TRAITS_FOR_CORRELATION <- 30 # Number of top phenotypes to include in correlation analysis

# --- CONFIGURABLE COLUMN NAMES IN YOUR INPUT TSV.GZ FILES ---
# !!! IMPORTANT: Adjust these to match the actual column names in your files !!!
PHENOTYPE_COL_NAME <- "match_file"  
SCORE_COL_NAME <- "product_total"  

# --- HELPER FUNCTION TO EXTRACT INDIVIDUAL ID FROM FILENAME ---
# This function extracts the ID used for SampleID (e.g., D1T0, D6T0a, D6T0b)
extract_individual_id <- function(filepath) {
  filename <- basename(filepath)
  # Split the filename by underscore and take the second element.
  id <- str_split_i(filename, "_", 2)
  return(id)
}

# --- HELPER FUNCTION TO CLEAN TRAIT LABELS (Python-style) ---
clean_label_python <- function(label) {
  label_str <- as.character(label)
  label_str <- str_replace(label_str, "_filtered\\.tsv\\.bgz", "")
  parts <- str_split(label_str, "-both_sexes", simplify = TRUE)
  return(parts[, 1])
}

# --- NEW HELPER FUNCTION: Crop string to the Nth underscore ---
crop_to_nth_underscore <- function(text, n) {
  if (is.na(text)) return(NA)
  parts <- str_split(text, "_")[[1]]
  if (length(parts) <= n) {
    return(text)
  } else {
    return(paste(parts[1:n], collapse = "_"))
  }
}

# --- MAIN DATA PROCESSING FUNCTION FOR A SINGLE RAW SCORE FILE ---
# Reads one .tsv.gz file, cleans phenotype labels, and attaches SampleID
process_single_raw_score_file <- function(filepath, phenotype_col_name, score_col_name, phenocodes_df, verbose = TRUE) {
  
  sample_id <- extract_individual_id(filepath)
  
  tryCatch({
    df_raw <- read_tsv(filepath, show_col_types = FALSE)
    
    required_cols <- c(phenotype_col_name, score_col_name)
    if (!all(required_cols %in% colnames(df_raw))) {
      stop(paste0("Missing required columns (", paste(required_cols, collapse = ", "), ") in ", filepath,
                  "\n  Found columns: ", paste(colnames(df_raw), collapse = ", ")))
    }
    
    df_processed <- df_raw %>%
      select(match_file = !!sym(phenotype_col_name), Score = !!sym(score_col_name)) %>%
      filter(!is.na(Score) & !is.infinite(Score) & !is.nan(Score))
    
    if (nrow(df_processed) == 0) {
      if (verbose) message(paste("  - Warning: No valid score data in '", filepath, "' after cleaning. Skipping."))
      return(NULL)
    }
    
    # Generate the base cleaned label (e.g., "UKB_B_N3_FSH")
    df_processed <- df_processed %>%
      mutate(
        label_clean_base = clean_label_python(match_file),
        # Apply the new cropping rule to the base label for the code part of the final Phenotype_Label
        cropped_code_label = sapply(label_clean_base, crop_to_nth_underscore, n = 3),
        filename_clean = str_replace(match_file, "_filtered\\.tsv\\.bgz", ".tsv.bgz")
      )
    
    # Merge with phenocodes for descriptive labels
    merged_df <- df_processed %>%
      left_join(phenocodes_df, by = c("filename_clean" = "filename")) %>%
      mutate(
        # Construct the final Phenotype_Label using the cropped code and description
        Phenotype_Label = ifelse(
          !is.na(description),
          paste0(cropped_code_label, " → ", description), # Cropped code + description
          cropped_code_label  # Fallback to cropped code if no description
        )
      )
    
    # Handle duplicate trait labels by taking mean score (if any)
    if (any(duplicated(merged_df$Phenotype_Label))) {
      if (verbose) message(paste0("  - Warning: Duplicate trait labels found in ", filepath, ". Taking mean score for duplicates."))
      merged_df <- merged_df %>%
        group_by(Phenotype_Label) %>%
        summarise(Score = mean(Score, na.rm = TRUE), .groups = 'drop')
    }
    
    merged_df$SampleID <- sample_id
    
    return(merged_df %>% select(Phenotype_Label, Score, SampleID))
    
  }, error = function(e) {
    if (verbose) message(paste("  - An error occurred while reading file", filepath, ":", e$message))
    return(NULL)
  })
}

# --- LOAD PHENOCODES ---
message(paste0("\n--- Loading phenocodes from: ", phenocodes_path, " ---"))
tryCatch({
  phenos <- read_tsv(phenocodes_path, show_col_types = FALSE)
  message(paste0("  - Successfully loaded phenocodes. Rows: ", nrow(phenos), ", Cols: ", ncol(phenos)))
  if (!all(c('filename', 'description') %in% colnames(phenos))) {
    stop(paste0("Error: 'filename' or 'description' column missing in '", phenocodes_path, "'. Please check phenocodes file."))
  }
}, error = function(e) {
  stop(paste0("Error loading phenocodes: ", e$message))
})

# --- MAIN PROCESSING: Read All Data ---
message("\n--- Reading all raw score data from input files ---")
all_raw_data_long <- data.frame()

if (length(input_files) == 0) {
  stop("Error: No input files specified in the 'input_files' list. Exiting.")
}
if (length(input_files) %% 2 != 0) {
  stop("Error: Input files list must contain an even number of files for paired analysis (C, O, C, O...).")
}

for (filepath in input_files) {
  if (file.exists(filepath)) {
    message(paste0("  - Processing '", filepath, "'..."))
    processed_df <- process_single_raw_score_file(
      filepath = filepath,
      phenotype_col_name = PHENOTYPE_COL_NAME,
      score_col_name = SCORE_COL_NAME,
      phenocodes_df = phenos,
      verbose = TRUE  
    )
    if (!is.null(processed_df)) {
      all_raw_data_long <- bind_rows(all_raw_data_long, processed_df)
    }
  } else {
    message(paste0("Warning: Input file not found: '", filepath, "'. Skipping."))
  }
}

if (nrow(all_raw_data_long) == 0) {
  stop("No data loaded from any input files. Please check file paths and ensure required columns are present and correctly named.")
}

# --- Filter out specific phenotypes ---
message("\n--- Filtering out unwanted phenotypes ---")
initial_pheno_count <- length(unique(all_raw_data_long$Phenotype_Label))

# Define the pattern for terms to remove, separated by '|' for OR logic in regex
terms_to_remove_pattern <- "Place of birth|television|Nap|Getting up|Z34|Alcohol|Time spent using computer|Hot drink|Lifetime number of sexual partners|fruit|Lamb|Worry|Cereal|PM: final answer"

all_raw_data_long <- all_raw_data_long %>%
  filter(!str_detect(Phenotype_Label, terms_to_remove_pattern))

filtered_pheno_count <- length(unique(all_raw_data_long$Phenotype_Label))
message(paste0("  - Removed phenotypes containing: '", gsub("\\|", "', '", terms_to_remove_pattern), "'. Initial unique phenotypes: ", initial_pheno_count, ", After filter: ", filtered_pheno_count, "."))


message(paste0("  - Total raw data rows loaded: ", nrow(all_raw_data_long)))

# --- Dumbbell Plot Generation ---
message("\n--- Generating individual dumbbell plots for each donor pair ---")
plot_list <- list()

for (i in seq(1, length(input_files), by = 2)) {
  file_c_path <- input_files[i]
  file_o_path <- input_files[i+1]
  
  sample_c_id <- extract_individual_id(file_c_path) # e.g., D1T0, D6T0a
  sample_o_id <- extract_individual_id(file_o_path) # e.g., D1T4, D6T4a
  
  # >>> ROBUST DONOR_ID EXTRACTION FOR ALL PLOTS/ANALYSES <<<
  d_part <- str_extract(sample_c_id, "D\\d+") 
  ab_part <- str_extract(sample_c_id, "[ab]") 
  
  if (!is.na(ab_part)) {
    donor_id <- paste0(d_part, ab_part) 
  } else {
    donor_id <- d_part
  }
  # --- END ROBUST DONOR_ID EXTRACTION ---
  
  message(paste0("  - Debugging: Processing pair: C_ID='", sample_c_id, "', O_ID='", sample_o_id, "'. Derived donor_id: '", donor_id, "'"))
  
  # Filter data for the current pair
  current_pair_data <- all_raw_data_long %>%
    filter(SampleID %in% c(sample_c_id, sample_o_id)) %>%
    pivot_wider(names_from = SampleID, values_from = Score) %>%
    filter(!is.na(.data[[sample_c_id]]) & !is.na(.data[[sample_o_id]])) %>%
    mutate(
      Score_C = .data[[sample_c_id]],
      Score_O = .data[[sample_o_id]],
      absolute_difference = abs(Score_O - Score_C)
    )
  
  if (nrow(current_pair_data) == 0) {
    message(paste0("    - No complete data for donor pair ", donor_id, ". Skipping plot."))
    plot_list[[paste0("plot_", donor_id)]] <- ggplot() + 
      labs(title = paste0(donor_id, " (No Data)"), 
           x = "Raw PGS Score", y = "Phenotype") + 
      theme_void() + 
      theme(plot.title = element_text(hjust = 0.5, size = 12, face = "italic"))
    next
  }
  
  # --- Per-Pair Trait Selection for Dumbbell Plots ---
  message(paste0("    - Selecting top phenotypes for pair ", donor_id, " (", PERCENTILE_THRESHOLD, "th percentile for dumbbell plot)..."))
  
  percentile_value_per_pair <- quantile(current_pair_data$absolute_difference, 1 - (PERCENTILE_THRESHOLD / 100.0), na.rm = TRUE)
  
  selected_top_phenotypes_for_this_pair <- current_pair_data %>%
    filter(absolute_difference >= percentile_value_per_pair) %>%
    arrange(desc(absolute_difference)) %>% 
    head(MAX_TRAITS_TO_PLOT_DUMBBELL) %>%
    pull(Phenotype_Label)
  
  if (length(selected_top_phenotypes_for_this_pair) == 0) {
    message(paste0("    - No phenotypes met the ", PERCENTILE_THRESHOLD, "th percentile threshold for pair ", donor_id, ". Selecting top 1 trait as fallback."))
    selected_top_phenotypes_for_this_pair <- head(current_pair_data %>% arrange(desc(absolute_difference)) %>% pull(Phenotype_Label), 1)
  }
  
  message(paste0("    - Selected ", length(selected_top_phenotypes_for_this_pair), " top phenotypes for pair ", donor_id, "."))
  
  pair_plot_data_filtered <- current_pair_data %>%
    filter(Phenotype_Label %in% selected_top_phenotypes_for_this_pair)
  
  if (nrow(pair_plot_data_filtered) == 0) {
    message(paste0("    - No data remaining after per-pair trait selection for donor pair ", donor_id, ". Skipping plot."))
    plot_list[[paste0("plot_", donor_id)]] <- ggplot() + 
      labs(title = paste0(donor_id, " (No Data)"), 
           x = "Raw PGS Score", y = "Phenotype") + 
      theme_void() + 
      theme(plot.title = element_text(hjust = 0.5, size = 12, face = "italic"))
    next
  }
  
  phenotype_order_for_plot <- pair_plot_data_filtered %>%
    arrange(Score_C) %>%
    pull(Phenotype_Label)
  
  pair_plot_data_filtered$Phenotype_Label <- factor(pair_plot_data_filtered$Phenotype_Label, levels = phenotype_order_for_plot)
  
  plot_points_data <- pair_plot_data_filtered %>%
    pivot_longer(
      cols = c(Score_C, Score_O),
      names_to = "Score_Type",
      values_to = "Score_Value"
    ) %>%
    mutate(
      Dot_Color_Category = case_when( 
        Score_Type == "Score_C" ~ "Control",
        Score_Type == "Score_O" ~ "Outcome"
      )
    )
  
  plot_min_score <- min(plot_points_data$Score_Value, na.rm = TRUE)
  plot_max_score <- max(plot_points_data$Score_Value, na.rm = TRUE)
  
  padding <- (plot_max_score - plot_min_score) * 0.1
  if (padding == 0) padding <- 1 
  plot_min_score_padded <- plot_min_score - padding
  plot_max_score_padded <- plot_max_score + padding
  
  if (plot_min_score_padded == plot_max_score_padded) {
    plot_min_score_padded <- plot_min_score - 0.1
    plot_max_score_padded <- plot_max_score + 0.1
  }
  
  p <- ggplot(pair_plot_data_filtered, aes(y = Phenotype_Label)) +
    geom_segment(aes(x = Score_C, xend = Score_O), 
                 linewidth = 0.8, color = "grey50") +
    geom_point(data = plot_points_data, 
               aes(x = Score_Value, color = Dot_Color_Category, shape = Dot_Color_Category), 
               size = 3, alpha=0.6) + 
    labs(
      subtitle = paste0("Donor: ", donor_id), 
      x = "Raw PGS Score",
      y = "Phenotype"
    ) +
    theme_minimal() +
    theme(
      plot.subtitle = element_text(hjust = 0.5, size = 12, face = "bold"),
      axis.title.x = element_text(size = 10, face = "bold", color = "black"),
      axis.title.y = element_text(size = 10, face = "bold", color = "black"),
      axis.text.x = element_text(size = 8, color = "black"),
      axis.text.y = element_text(size = 7, color = "black"),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.title = element_text(size = 9, face = "bold"),
      legend.text = element_text(size = 8),
      panel.grid.major.y = element_line(color = "lightgray", linetype = "dotted", linewidth = 0.5),
      panel.grid.major.x = element_line(color = "lightgray", linetype = "dotted", linewidth = 0.5),
      panel.grid.minor = element_blank(),
      plot.margin = unit(c(0.25, 0.25, 0.25, 0.25), "cm")
    ) +
    coord_cartesian(xlim = c(plot_min_score_padded, plot_max_score_padded)) +
    scale_color_manual(
      name = "Sample Type",
      values = c("Control" = "red", "Outcome" = "blue"),
      labels = c(paste0("Control (", sample_c_id, ")"), paste0("Outcome (", sample_o_id, ")"))
    ) +
    scale_shape_manual(
      name = "Sample Type", 
      values = c("Control" = 16, "Outcome" = 17), 
      labels = c(paste0("Control (", sample_c_id, ")"), paste0("Outcome (", sample_o_id, ")")) 
    )
  
  list_key <- paste0("plot_", donor_id)
  message(paste0("  - Debugging: Storing plot for donor_id '", donor_id, "' with key '", list_key, "'"))
  plot_list[[list_key]] <- p
}

message("\n--- Debugging: Keys in plot_list after generation (Dumbbell Plots) ---")
print(names(plot_list))
message(paste0("  - Total dumbbell plots generated (by unique keys): ", length(plot_list)))

# --- Create Composite Dumbbell Plot ---
message("\n--- Creating composite dumbbell plot ---")

n_plots_dumbbell <- length(plot_list)
n_cols_composite_dumbbell <- 2 
n_rows_composite_dumbbell <- ceiling(n_plots_dumbbell / n_cols_composite_dumbbell) 

composite_dumbbell_plot <- wrap_plots(plot_list, ncol = n_cols_composite_dumbbell) +
  plot_annotation(
    title = paste0("PGS Raw Score Changes by Donor\n(Top ", MAX_TRAITS_TO_PLOT_DUMBBELL, " Phenotypes per Pair, displaying 99th percentile)"),
    theme = theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold", margin = margin(b = 10))) 
  ) & theme(legend.position = "bottom") 

height_per_panel_dumbbell <- max(7, 0.25 * MAX_TRAITS_TO_PLOT_DUMBBELL) 
composite_dumbbell_plot_height_in <- n_rows_composite_dumbbell * height_per_panel_dumbbell
composite_dumbbell_plot_width_in <- n_cols_composite_dumbbell * 8 

if (composite_dumbbell_plot_height_in > 40) composite_dumbbell_plot_height_in <- 40  
if (composite_dumbbell_plot_width_in > 30) composite_dumbbell_plot_width_in <- 30  

ggsave(output_composite_dumbbell_plot_filename,
       plot = composite_dumbbell_plot,
       width = composite_dumbbell_plot_width_in,
       height = composite_dumbbell_plot_height_in,
       units = "in", dpi = 300, bg = "white")

message(paste0("\n--- Composite dumbbell plot saved to '", output_composite_dumbbell_plot_filename, "' ---"))


# --- Cross-Phenotype Correlation Analysis of Change Scores ---
message("\n--- Starting Cross-Phenotype Correlation Analysis of Change Scores ---")

# 1. Prepare data for calculating change scores for each donor and phenotype
message("  - Restructuring data to calculate change scores per donor and phenotype...")
df_paired_scores <- all_raw_data_long %>%
  mutate(
    d_part = str_extract(SampleID, "D\\d+"),
    ab_part = str_extract(SampleID, "[ab]"),
    donor_id = ifelse(!is.na(ab_part), paste0(d_part, ab_part), d_part)
  ) %>%
  mutate(
    Timepoint = case_when(
      str_detect(SampleID, "T0") ~ "T0",
      str_detect(SampleID, "T4") ~ "T4",
      TRUE ~ NA_character_
    )
  ) %>%
  pivot_wider(
    id_cols = c(donor_id, Phenotype_Label),
    names_from = Timepoint,
    values_from = Score,
    names_prefix = "Score_"
  ) %>%
  filter(!is.na(Score_T0) & !is.na(Score_T4)) %>%
  mutate(Change_Score = Score_T4 - Score_T0)

if (nrow(df_paired_scores) == 0) {
  stop("Error: No complete paired T0 and T4 data found for correlation analysis after restructuring.")
}
message(paste0("  - Data restructured. Total paired change score entries: ", nrow(df_paired_scores)))

# 2. Select Global Top Phenotypes for Correlation
message("  - Selecting top phenotypes for global correlation analysis...")
global_phenotype_avg_change <- df_paired_scores %>%
  group_by(Phenotype_Label) %>%
  summarise(
    Average_Abs_Change = mean(abs(Change_Score), na.rm = TRUE),
    N_Donors = n(),
    .groups = 'drop'
  ) %>%
  arrange(desc(Average_Abs_Change))

selected_phenotypes_for_correlation <- global_phenotype_avg_change %>%
  filter(N_Donors >= 2) %>% 
  head(MAX_TRAITS_FOR_CORRELATION) %>%
  pull(Phenotype_Label)

if (length(selected_phenotypes_for_correlation) < 2) {
  stop(paste0("Error: Not enough phenotypes (need at least 2) with sufficient data for correlation analysis. Found ", length(selected_phenotypes_for_correlation), " selected."))
}
message(paste0("  - Selected ", length(selected_phenotypes_for_correlation), " top phenotypes for correlation."))

# 3. Pivot data wider for correlation matrix (rows = donors, columns = phenotype change scores)
df_change_scores_wide <- df_paired_scores %>%
  filter(Phenotype_Label %in% selected_phenotypes_for_correlation) %>%
  select(donor_id, Phenotype_Label, Change_Score) %>%
  pivot_wider(
    names_from = Phenotype_Label,
    values_from = Change_Score
  ) %>%
  column_to_rownames("donor_id") 

message(paste0("  - Prepared correlation matrix data for ", nrow(df_change_scores_wide), " donors and ", ncol(df_change_scores_wide), " phenotypes."))

# 4. Perform Correlation Analysis using Hmisc::rcorr
message("  - Calculating Pearson correlation matrix and p-values...")
correlation_results <- rcorr(as.matrix(df_change_scores_wide), type = "pearson")

R_matrix <- correlation_results$r
P_matrix <- correlation_results$P

# 5. Report Significant Correlations (p < 0.05)
message("\n--- Significant Cross-Phenotype Correlations of Change Scores (p < 0.05) ---")
significant_correlations <- data.frame(
  Phenotype1 = rownames(R_matrix)[row(R_matrix)],
  Phenotype2 = colnames(R_matrix)[col(R_matrix)],
  Pearson_r = as.vector(R_matrix),
  P_value = as.vector(P_matrix)
) %>%
  filter(Phenotype1 != Phenotype2) %>% 
  rowwise() %>%
  mutate(
    Pair_ID = paste(sort(c(Phenotype1, Phenotype2)), collapse = "---")
  ) %>%
  ungroup() %>%
  distinct(Pair_ID, .keep_all = TRUE) %>% 
  select(-Pair_ID) %>%
  filter(P_value < 0.05) %>% 
  arrange(P_value) 

if (nrow(significant_correlations) > 0) {
  print(significant_correlations)
  message("\nNote: A small number of donors (N=", nrow(df_change_scores_wide), ") may lead to unstable correlations and high p-values. Interpret with caution.")
} else {
  message("  - No significant correlations (p < 0.05) found among the selected phenotypes.")
}

# 6. Visualize Correlation Matrix (Heatmap)
message("\n--- Generating Correlation Heatmap ---")

heatmap_title <- paste0("Pearson Correlation of PGS Score Changes Across Donors (N=", nrow(df_change_scores_wide), ")\nTop ", length(selected_phenotypes_for_correlation), " Phenotypes by Avg. Absolute Change")

png(output_correlation_heatmap_filename, 
    width = 10 + ncol(R_matrix)*0.2, 
    height = 10 + ncol(R_matrix)*0.2, 
    units = "in", res = 300)

corrplot(R_matrix, 
         p.mat = P_matrix, 
         sig.level = 0.05, 
         insig = "blank", 
         method = "color", 
         type = "upper", 
         tl.col = "black", 
         tl.srt = 45, 
         tl.cex = 0.6, 
         number.cex = 0.5, 
         addCoef.col = "black", 
         col = colorRampPalette(c("blue", "white", "red"))(200), 
         title = heatmap_title,
         mar = c(0,0,3,0)
)

dev.off() 
message(paste0("--- Correlation heatmap saved to '", output_correlation_heatmap_filename, "' ---"))

# 7. Generate a Single Scatter Plot for the Most Negatively Correlated Pair
message("\n--- Generating a single scatter plot for the most negatively correlated change scores ---")

# Prioritize significant negative correlations, then fall back to most negative overall
most_negative_correlation <- significant_correlations %>%
  filter(Pearson_r < 0) %>%
  arrange(Pearson_r) %>%
  head(1)

if (nrow(most_negative_correlation) == 0) {
  # If no *significant* negative correlations, find the most negative one regardless of significance
  message("  - No *significant* negative correlations found. Looking for the overall most negative correlation.")
  most_negative_correlation <- data.frame(
    Phenotype1 = rownames(R_matrix)[row(R_matrix)],
    Phenotype2 = colnames(R_matrix)[col(R_matrix)],
    Pearson_r = as.vector(R_matrix),
    P_value = as.vector(P_matrix)
  ) %>%
    filter(Phenotype1 != Phenotype2) %>%
    rowwise() %>%
    mutate(Pair_ID = paste(sort(c(Phenotype1, Phenotype2)), collapse = "---")) %>%
    ungroup() %>%
    distinct(Pair_ID, .keep_all = TRUE) %>%
    select(-Pair_ID) %>%
    arrange(Pearson_r) %>% # Sort by most negative R
    head(1)
}

if (nrow(most_negative_correlation) > 0) {
  pheno1_label <- most_negative_correlation$Phenotype1[1]
  pheno2_label <- most_negative_correlation$Phenotype2[1]
  pearson_r <- most_negative_correlation$Pearson_r[1]
  p_value <- most_negative_correlation$P_value[1]
  
  plot_data_scatter <- df_change_scores_wide %>%
    select(all_of(c(pheno1_label, pheno2_label))) %>%
    rename(Pheno1_Change = all_of(pheno1_label),
           Pheno2_Change = all_of(pheno2_label)) %>%
    mutate(donor_id_col = rownames(.)) 
  
  p_scatter_single <- ggplot(plot_data_scatter, aes(x = Pheno1_Change, y = Pheno2_Change)) +
    geom_point(size = 4, color = "darkgreen") +
    geom_text(aes(label = donor_id_col), hjust = -0.3, vjust = 0.5, size = 3.5, color = "darkred") + 
    geom_smooth(method = "lm", se = FALSE, color = "blue", linetype = "dashed") + 
    labs(
      title = paste0("Change Score Correlation: ", pheno1_label, " vs ", pheno2_label),
      subtitle = paste0("Pearson r = ", round(pearson_r, 3), ", p = ", format.pval(p_value, digits = 3),
                        "\nEach point represents a donor. Note: Small N (", nrow(df_change_scores_wide), ") makes correlations highly variable."),
      x = paste0("Change in ", pheno1_label),
      y = paste0("Change in ", pheno2_label)
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 9, face = "italic", color = "darkgray"),
      axis.title = element_text(size = 10, face = "bold"),
      axis.text = element_text(size = 8)
    )
  
  ggsave(output_single_scatter_plot_filename,
         plot = p_scatter_single,
         width = 8, height = 6,
         units = "in", dpi = 300, bg = "white")
  
  message(paste0("--- Single scatter plot saved to '", output_single_scatter_plot_filename, "' ---"))
  
} else {
  message("  - No negative correlations found at all to plot a scatter plot.")
}

message("\n--- Analysis complete ---")



# --- EXPORT ALL CORRELATIONS TO CSV ---
message("\n--- Exporting all Pearson correlations and p-values to CSV ---")

# Combine R_matrix and P_matrix into a single data frame
# First, convert matrices to long format
all_correlations_r_long <- as.data.frame(R_matrix) %>%
  tibble::rownames_to_column(var = "Phenotype1") %>%
  pivot_longer(cols = -Phenotype1, names_to = "Phenotype2", values_to = "Pearson_r")

all_correlations_p_long <- as.data.frame(P_matrix) %>%
  tibble::rownames_to_column(var = "Phenotype1") %>%
  pivot_longer(cols = -Phenotype1, names_to = "Phenotype2", values_to = "P_value")

# Merge the R and P values
all_correlations_table_export <- left_join(all_correlations_r_long, all_correlations_p_long, 
                                           by = c("Phenotype1", "Phenotype2"))

# Filter out self-correlations (where Phenotype1 == Phenotype2)
all_correlations_table_export <- all_correlations_table_export %>%
  filter(Phenotype1 != Phenotype2)

# Ensure unique pairs (e.g., A-B is the same as B-A, keep only one entry)
all_correlations_table_export <- all_correlations_table_export %>%
  rowwise() %>%
  mutate(Pair_ID = paste(sort(c(Phenotype1, Phenotype2)), collapse = "---")) %>%
  ungroup() %>%
  distinct(Pair_ID, .keep_all = TRUE) %>%
  select(-Pair_ID) %>% # Remove the temporary Pair_ID column
  arrange(P_value) # Optional: Sort by P_value for easier review

# Define the output filename
output_all_correlations_filename <- "all_phenotype_change_correlations_R_SU.csv"

# Write the data frame to a CSV file
write.csv(all_correlations_table_export, file = output_all_correlations_filename, row.names = FALSE)

message(paste0("--- All correlation data (Pearson r and P values) saved to '", output_all_correlations_filename, "' ---"))