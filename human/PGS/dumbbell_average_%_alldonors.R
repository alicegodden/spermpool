Title: Absolute raw average, and plot of percentage change across all donors 
Author: Dr. Alice M. Godden

# SwimUp
# ----------------------------------------------------------------------
# PGS Dumbbell Plots: Percentage Change (Delta %) FACETED by Delta % Magnitude
# FIX: Y-axis (Phenotype Labels) is now independent for each facet, 
# ensuring only labels with data points are displayed.
# ----------------------------------------------------------------------

# --- LIBRARIES ---
library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(ggplot2)
library(tools)
library(scales)   # For percent_format

# --- CONFIGURATION PARAMETERS ---
input_files <- c(
  "GRCh37_D1T0_merged.dedup_bcftools__CHR_nc_to_chr.nochr.aaf_scores.tsv.gz",
  "GRCh37_D1T4_merged.dedup_RG_bcftools__CHR_nc_to_chr.nochr.aaf_scores.tsv.gz",
  "GRCh37_D2T0_merged.dedup_bcftools__CHR_nc_to_chr.nochr.aaf_scores.tsv.gz",
  "GRCh37_D2T4_merged.dedup_bcftools__CHR_nc_to_chr.nochr.aaf_scores.tsv.gz"
)
phenocodes_path <- "phenocodes"  # must contain columns: filename, description
PHENOTYPE_COL_NAME <- "match_file"
SCORE_COL_NAME <- "product_total"
PERCENTILE <- 99
TOP_N <- 50
DELTA_PERCENT_BANDS <- 4 # Number of facets based on Delta % magnitude (quartiles)

# --- HELPER FUNCTIONS ---

extract_individual_id <- function(filepath) {
  filename <- basename(filepath)
  str_extract(filename, "D\\d+T\\d+[ab]?")
}

clean_label_python <- function(label) {
  label_str <- as.character(label)
  label_str <- str_replace(label_str, "_filtered\\.tsv\\.bgz", "")
  parts <- str_split(label_str, "-both_sexes", simplify = TRUE)
  return(parts[, 1])
}

crop_to_nth_underscore <- function(text, n) {
  if (is.na(text)) return(NA)
  parts <- str_split(text, "_")[[1]]
  if (length(parts) <= n) {
    return(text)
  } else {
    return(paste(parts[1:n], collapse = "_"))
  }
}

process_single_raw_score_file <- function(filepath, phenotype_col_name, score_col_name, phenocodes_df) {
  sample_id <- extract_individual_id(filepath)
  
  df_raw <- read_tsv(filepath, show_col_types = FALSE) %>%
    select(match_file = !!sym(phenotype_col_name), Score = !!sym(score_col_name)) %>%
    filter(!is.na(Score) & !is.infinite(Score) & !is.nan(Score))
  
  df_processed <- df_raw %>%
    mutate(
      label_clean_base = clean_label_python(match_file),
      cropped_code_label = sapply(label_clean_base, crop_to_nth_underscore, n = 3),
      filename_clean = str_replace(match_file, "_filtered\\.tsv\\.bgz", ".tsv.bgz")
    )
  
  merged_df <- df_processed %>%
    left_join(phenocodes_df, by = c("filename_clean" = "filename")) %>%
    mutate(
      Phenotype_Label = ifelse(!is.na(description), description, cropped_code_label)
    )
  
  merged_df <- merged_df %>%
    group_by(Phenotype_Label) %>%
    summarise(Score = mean(Score, na.rm = TRUE), .groups = 'drop')
  
  merged_df$SampleID <- sample_id
  
  return(merged_df %>% select(Phenotype_Label, Score, SampleID))
}


# --- PLOTTING FUNCTION (Adapted for Percentage Change and Faceting) ---
make_dumbbell_plot_delta <- function(df, plot_title="") {
  
  plot_points_data <- df %>%
    pivot_longer(cols = c(Score_C_perc, Score_O_perc),
                 names_to = "Score_Type",
                 values_to = "Score_Value") %>%
    mutate(Dot_Color = ifelse(Score_Type=="Score_C_perc", "T0","T4"))
  
  p <- ggplot(df, aes(y=Phenotype_Label)) +
    geom_segment(aes(x=Score_C_perc, xend=Score_O_perc),
                 color="grey50", linewidth=0.8) +
    geom_point(data=plot_points_data,
               aes(x=Score_Value, color=Dot_Color, shape=Dot_Color),
               size=3) +
    # --- KEY FIX: Use facet_wrap with scales="free" and ncol=1 ---
    # scales="free" makes both X (score) and Y (label) axes independent.
    # ncol=1 stacks the facets vertically.
    facet_wrap(~Score_Group, ncol=1, scales="free") + 
    labs(title=plot_title,
         x="Percentage Change from T0 (%)",
         y="Phenotype") +
    scale_color_manual(values=c("T0"="red","T4"="blue"), name="Time Point") +
    scale_shape_manual(values=c("T0"=16,"T4"=17), name="Time Point") +
    # Format X-axis as percentage and add a vertical line at 0%
    scale_x_continuous(labels = scales::percent_format(scale=1)) +
    geom_vline(xintercept=0, linetype="dashed", color="darkgrey", linewidth=0.5) +
    theme_minimal() +
    theme(
      plot.title=element_text(hjust=0.5, face="bold"),
      strip.text.x = element_text(size=9, face="bold"),
      strip.background = element_rect(fill="grey95", color="grey50"),
      axis.text.y=element_text(size=9),
      legend.position = "bottom"
    )
  
  return(p)
}


# --- MAIN SCRIPT EXECUTION ---

# --- LOAD PHENOCODES ---
message(paste0("\n--- Loading phenocodes from: ", phenocodes_path, " ---"))
tryCatch({
  phenos <- read_tsv(phenocodes_path, show_col_types = FALSE)
  if (!all(c('filename', 'description') %in% colnames(phenos))) {
    stop(paste0("Error: 'filename' or 'description' column missing in '", phenocodes_path, "'."))
  }
}, error = function(e) {
  stop(paste0("Error loading phenocodes: ", e$message))
})

# --- LOAD ALL DATA ---
all_data <- data.frame()
for (f in input_files) {
  if(file.exists(f)) {
    df <- process_single_raw_score_file(
      filepath = f,
      phenotype_col_name = PHENOTYPE_COL_NAME,
      score_col_name = SCORE_COL_NAME,
      phenocodes_df = phenos
    )
    all_data <- bind_rows(all_data, df)
  } else {
    message(paste("File not found:", f))
  }
}

if (nrow(all_data) == 0) {
  stop("No data loaded. Check file paths and phenocodes file.")
}

# --- GLOBAL CHANGES (for filtering) ---
global_changes <- data.frame()

for (i in seq(1,length(input_files), by=2)) {
  file_c <- input_files[i]
  file_o <- input_files[i+1]
  
  id_c <- extract_individual_id(file_c)
  id_o <- extract_individual_id(file_o)
  
  df_pair <- all_data %>%
    filter(SampleID %in% c(id_c,id_o)) %>%
    pivot_wider(names_from = SampleID,
                values_from = Score,
                values_fn = list(Score = mean),
                values_fill = NA) %>%
    filter(!is.na(.data[[id_c]]) & !is.na(.data[[id_o]])) %>%
    mutate(
      Score_C = .data[[id_c]],
      Score_O = .data[[id_o]],
      abs_diff = abs(Score_O - Score_C)
    )
  
  global_changes <- bind_rows(global_changes,
                              df_pair %>% select(Phenotype_Label, abs_diff))
}

# --- GLOBAL TOP PHENOTYPES (Filtering) ---
threshold <- quantile(global_changes$abs_diff, PERCENTILE/100, na.rm=TRUE)
top_traits <- global_changes %>%
  filter(abs_diff >= threshold) %>%
  distinct(Phenotype_Label) %>%
  head(TOP_N) %>%
  pull(Phenotype_Label)

# --- GLOBAL AVERAGE SCORES (Delta % Calculation and Grouping) ---
avg_df <- all_data %>%
  filter(Phenotype_Label %in% top_traits) %>%
  mutate(Time = ifelse(str_detect(SampleID, "T0"), "C","O")) %>%
  group_by(Phenotype_Label, Time) %>%
  summarise(mean_score = mean(Score, na.rm=TRUE), .groups="drop") %>%
  pivot_wider(names_from = Time, values_from = mean_score) %>%
  
  # Filter out rows where C or O are NA or C is 0 (ensures valid percentage calculation)
  filter(!is.na(C) & !is.na(O) & C != 0) %>%
  
  mutate(
    Score_C_perc = 0,
    Score_O_perc = ((O - C) / C) * 100 
  ) %>%
  
  # --- GROUPING: Facet by Delta % Magnitude (Quartiles) ---
  mutate(Score_Group_Temp = ntile(Score_O_perc, DELTA_PERCENT_BANDS)) %>%
  group_by(Score_Group_Temp) %>%  
  mutate(
    Score_Group = paste0("Quartile ", Score_Group_Temp, " (", 
                         format(min(Score_O_perc), digits = 2), 
                         "% to ", 
                         format(max(Score_O_perc), digits = 2), "%)")
  ) %>%
  ungroup() %>%
  
  # --- Final ordering and factoring for plotting ---
  # 1. Order the Score_Group factor correctly
  mutate(Score_Group = factor(Score_Group, 
                              levels = unique(Score_Group[order(Score_Group_Temp)]))) %>%
  group_by(Score_Group) %>%
  # 2. Sort Y-axis labels within each facet by the percentage change
  arrange(Score_Group, Score_O_perc) %>%
  
  # 3. Factor Y-axis only by the remaining unique Phenotype_Labels
  mutate(Phenotype_Label = factor(Phenotype_Label, levels = unique(Phenotype_Label))) %>%
  
  ungroup() %>%
  select(-Score_Group_Temp)


write.csv(avg_df %>% select(Phenotype_Label, C, O, Score_O_perc, Score_Group), 
          "PGS_GlobalTop50Average_DeltaPercentChange_Faceted.csv", row.names = FALSE)

# --- PLOT (ONE GLOBAL FACETED PLOT) ---
avg_plot <- make_dumbbell_plot_delta(avg_df, paste("Average % change of absolute raw score across donors: SwimUp (n=", DELTA_PERCENT_BANDS, ")", sep=""))

# --- SAVE PLOTS ---
ggsave("PGS_GlobalTop50Average_DeltaPercentChange_Faceted.png", avg_plot, width=12, height=14, dpi=300)

message("Script completed. The plot now uses `facet_wrap` with `scales=\"free\"` to ensure the Y-axis (phenotype labels) is unique to the data present in each facet.")



# Methylcellulose

# ----------------------------------------------------------------------
# PGS Dumbbell Plots: Percentage Change (Delta %) FACETED by Delta % Magnitude
# FIX: Y-axis (Phenotype Labels) is now independent for each facet, 
# ensuring only labels with data points are displayed.
# ----------------------------------------------------------------------

# --- LIBRARIES ---
library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(ggplot2)
library(tools)
library(scales)   # For percent_format

# --- CONFIGURATION PARAMETERS ---
input_files <- c(
  "GRCh37_M8_C1.dedup_bcftools__CHR_nc_to_chr.nochr.aaf_scores.tsv.gz",
  "GRCh37_M8_O1.dedup_bcftools__CHR_nc_to_chr.nochr.aaf_scores.tsv.gz",
  "GRCh37_M11_C2._merged_dedup_bcftools__CHR_nc_to_chr.nochr.aaf_scores.tsv.gz",
  "GRCh37_M11_O1.dedup_bcftools__CHR_nc_to_chr.nochr.aaf_scores.tsv.gz",
  "GRCh37_M12_C1.dedup_bcftools__CHR_nc_to_chr.nochr.aaf_scores.tsv.gz",
  "GRCh37_M12_O1.dedup_bcftools__CHR_nc_to_chr.nochr.aaf_scores.tsv.gz"
  "GRCh37_M13_C1.dedup_bcftools__CHR_nc_to_chr.nochr.aaf_scores.tsv.gz",
  "GRCh37_M13_O1.dedup_bcftools__CHR_nc_to_chr.nochr.aaf_scores.tsv.gz"
)
phenocodes_path <- "phenocodes"  # must contain columns: filename, description
PHENOTYPE_COL_NAME <- "match_file"
SCORE_COL_NAME <- "product_total"
PERCENTILE <- 99
TOP_N <- 50
DELTA_PERCENT_BANDS <- 4 # Number of facets based on Delta % magnitude (quartiles)

# --- HELPER FUNCTIONS ---

extract_individual_id <- function(filepath) {
  filename <- basename(filepath)
  str_extract(filename, "D\\d+T\\d+[ab]?")
}

clean_label_python <- function(label) {
  label_str <- as.character(label)
  label_str <- str_replace(label_str, "_filtered\\.tsv\\.bgz", "")
  parts <- str_split(label_str, "-both_sexes", simplify = TRUE)
  return(parts[, 1])
}

crop_to_nth_underscore <- function(text, n) {
  if (is.na(text)) return(NA)
  parts <- str_split(text, "_")[[1]]
  if (length(parts) <= n) {
    return(text)
  } else {
    return(paste(parts[1:n], collapse = "_"))
  }
}

process_single_raw_score_file <- function(filepath, phenotype_col_name, score_col_name, phenocodes_df) {
  sample_id <- extract_individual_id(filepath)
  
  df_raw <- read_tsv(filepath, show_col_types = FALSE) %>%
    select(match_file = !!sym(phenotype_col_name), Score = !!sym(score_col_name)) %>%
    filter(!is.na(Score) & !is.infinite(Score) & !is.nan(Score))
  
  df_processed <- df_raw %>%
    mutate(
      label_clean_base = clean_label_python(match_file),
      cropped_code_label = sapply(label_clean_base, crop_to_nth_underscore, n = 3),
      filename_clean = str_replace(match_file, "_filtered\\.tsv\\.bgz", ".tsv.bgz")
    )
  
  merged_df <- df_processed %>%
    left_join(phenocodes_df, by = c("filename_clean" = "filename")) %>%
    mutate(
      Phenotype_Label = ifelse(!is.na(description), description, cropped_code_label)
    )
  
  merged_df <- merged_df %>%
    group_by(Phenotype_Label) %>%
    summarise(Score = mean(Score, na.rm = TRUE), .groups = 'drop')
  
  merged_df$SampleID <- sample_id
  
  return(merged_df %>% select(Phenotype_Label, Score, SampleID))
}


# --- PLOTTING FUNCTION (Adapted for Percentage Change and Faceting) ---
make_dumbbell_plot_delta <- function(df, plot_title="") {
  
  plot_points_data <- df %>%
    pivot_longer(cols = c(Score_C_perc, Score_O_perc),
                 names_to = "Score_Type",
                 values_to = "Score_Value") %>%
    mutate(Dot_Color = ifelse(Score_Type=="Score_C_perc", "T0","T4"))
  
  p <- ggplot(df, aes(y=Phenotype_Label)) +
    geom_segment(aes(x=Score_C_perc, xend=Score_O_perc),
                 color="grey50", linewidth=0.8) +
    geom_point(data=plot_points_data,
               aes(x=Score_Value, color=Dot_Color, shape=Dot_Color),
               size=3) +
    # --- KEY FIX: Use facet_wrap with scales="free" and ncol=1 ---
    # scales="free" makes both X (score) and Y (label) axes independent.
    # ncol=1 stacks the facets vertically.
    facet_wrap(~Score_Group, ncol=1, scales="free") + 
    labs(title=plot_title,
         x="Percentage Change from T0 (%)",
         y="Phenotype") +
    scale_color_manual(values=c("T0"="red","T4"="blue"), name="Time Point") +
    scale_shape_manual(values=c("T0"=16,"T4"=17), name="Time Point") +
    # Format X-axis as percentage and add a vertical line at 0%
    scale_x_continuous(labels = scales::percent_format(scale=1)) +
    geom_vline(xintercept=0, linetype="dashed", color="darkgrey", linewidth=0.5) +
    theme_minimal() +
    theme(
      plot.title=element_text(hjust=0.5, face="bold"),
      strip.text.x = element_text(size=9, face="bold"),
      strip.background = element_rect(fill="grey95", color="grey50"),
      axis.text.y=element_text(size=9),
      legend.position = "bottom"
    )
  
  return(p)
}


# --- MAIN SCRIPT EXECUTION ---

# --- LOAD PHENOCODES ---
message(paste0("\n--- Loading phenocodes from: ", phenocodes_path, " ---"))
tryCatch({
  phenos <- read_tsv(phenocodes_path, show_col_types = FALSE)
  if (!all(c('filename', 'description') %in% colnames(phenos))) {
    stop(paste0("Error: 'filename' or 'description' column missing in '", phenocodes_path, "'."))
  }
}, error = function(e) {
  stop(paste0("Error loading phenocodes: ", e$message))
})

# --- LOAD ALL DATA ---
all_data <- data.frame()
for (f in input_files) {
  if(file.exists(f)) {
    df <- process_single_raw_score_file(
      filepath = f,
      phenotype_col_name = PHENOTYPE_COL_NAME,
      score_col_name = SCORE_COL_NAME,
      phenocodes_df = phenos
    )
    all_data <- bind_rows(all_data, df)
  } else {
    message(paste("File not found:", f))
  }
}

if (nrow(all_data) == 0) {
  stop("No data loaded. Check file paths and phenocodes file.")
}

# --- GLOBAL CHANGES (for filtering) ---
global_changes <- data.frame()

for (i in seq(1,length(input_files), by=2)) {
  file_c <- input_files[i]
  file_o <- input_files[i+1]
  
  id_c <- extract_individual_id(file_c)
  id_o <- extract_individual_id(file_o)
  
  df_pair <- all_data %>%
    filter(SampleID %in% c(id_c,id_o)) %>%
    pivot_wider(names_from = SampleID,
                values_from = Score,
                values_fn = list(Score = mean),
                values_fill = NA) %>%
    filter(!is.na(.data[[id_c]]) & !is.na(.data[[id_o]])) %>%
    mutate(
      Score_C = .data[[id_c]],
      Score_O = .data[[id_o]],
      abs_diff = abs(Score_O - Score_C)
    )
  
  global_changes <- bind_rows(global_changes,
                              df_pair %>% select(Phenotype_Label, abs_diff))
}

# --- GLOBAL TOP PHENOTYPES (Filtering) ---
threshold <- quantile(global_changes$abs_diff, PERCENTILE/100, na.rm=TRUE)
top_traits <- global_changes %>%
  filter(abs_diff >= threshold) %>%
  distinct(Phenotype_Label) %>%
  head(TOP_N) %>%
  pull(Phenotype_Label)

# --- GLOBAL AVERAGE SCORES (Delta % Calculation and Grouping) ---
avg_df <- all_data %>%
  filter(Phenotype_Label %in% top_traits) %>%
  mutate(Time = ifelse(str_detect(SampleID, "T0"), "C","O")) %>%
  group_by(Phenotype_Label, Time) %>%
  summarise(mean_score = mean(Score, na.rm=TRUE), .groups="drop") %>%
  pivot_wider(names_from = Time, values_from = mean_score) %>%
  
  # Filter out rows where C or O are NA or C is 0 (ensures valid percentage calculation)
  filter(!is.na(C) & !is.na(O) & C != 0) %>%
  
  mutate(
    Score_C_perc = 0,
    Score_O_perc = ((O - C) / C) * 100 
  ) %>%
  
  # --- GROUPING: Facet by Delta % Magnitude (Quartiles) ---
  mutate(Score_Group_Temp = ntile(Score_O_perc, DELTA_PERCENT_BANDS)) %>%
  group_by(Score_Group_Temp) %>%  
  mutate(
    Score_Group = paste0("Quartile ", Score_Group_Temp, " (", 
                         format(min(Score_O_perc), digits = 2), 
                         "% to ", 
                         format(max(Score_O_perc), digits = 2), "%)")
  ) %>%
  ungroup() %>%
  
  # --- Final ordering and factoring for plotting ---
  # 1. Order the Score_Group factor correctly
  mutate(Score_Group = factor(Score_Group, 
                              levels = unique(Score_Group[order(Score_Group_Temp)]))) %>%
  group_by(Score_Group) %>%
  # 2. Sort Y-axis labels within each facet by the percentage change
  arrange(Score_Group, Score_O_perc) %>%
  
  # 3. Factor Y-axis only by the remaining unique Phenotype_Labels
  mutate(Phenotype_Label = factor(Phenotype_Label, levels = unique(Phenotype_Label))) %>%
  
  ungroup() %>%
  select(-Score_Group_Temp)


write.csv(avg_df %>% select(Phenotype_Label, C, O, Score_O_perc, Score_Group), 
          "PGS_GlobalTop50Average_DeltaPercentChange_Faceted_MC.csv", row.names = FALSE)

# --- PLOT (ONE GLOBAL FACETED PLOT) ---
avg_plot <- make_dumbbell_plot_delta(avg_df, paste("Average % change of absolute raw score across donors: MethylCellulose (n=", DELTA_PERCENT_BANDS, ")", sep=""))

# --- SAVE PLOTS ---
ggsave("PGS_GlobalTop50Average_DeltaPercentChange_Faceted_MC.png", avg_plot, width=12, height=14, dpi=300)

message("Script completed. The plot now uses `facet_wrap` with `scales=\"free\"` to ensure the Y-axis (phenotype labels) is unique to the data present in each facet.")
