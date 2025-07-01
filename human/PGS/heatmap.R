# Heatmap of PGS Score Differences with Phenocode Labels and Individual-wise Z-score Normalization

# --- SET WORKING DIRECTORY ---
# <<< IMPORTANT: Adjust this path to where your input files and phenocodes file are located >>>
setwd("~/PycharmProjects/polygenic/PGS_scores_test/SU/T0") 

# --- LOAD NECESSARY LIBRARIES ---
library(dplyr)      # For data manipulation
library(readr)      # For reading data (read_tsv)
library(ggplot2)    # For plotting
library(stringr)    # For string manipulation
library(tools)      # For file path manipulation (basename)
library(tidyr)      # For pivot_longer, pivot_wider
library(purrr)      # For reduce (to combine data frames efficiently)

# --- CONFIGURATION ---
# >>> IMPORTANT: Ensure these match your actual file structure and desired behavior <<<

# Input files (score difference files from the Python script's 'MC' list)
#input_files <- c(
#  "PGSscores_difference_MC_M8_C-O.txt",
#  "PGSscores_difference_MC_M11_C-O.txt",
#  "PGSscores_difference_MC_M12_C-O.txt",
#  "PGSscores_difference_MC_M13_C-O.txt"
#)

input_files <- c(
  "PGSscores_difference_D1T4-D1T0.txt",
  "PGSscores_difference_D2T4-D2T0.txt",
  "PGSscores_difference_D4T4-D4T0.txt",
  "PGSscores_difference_D6aT4-D6T0a.txt",
  "PGSscores_difference_D6bT4-D6bT0.txt"
)

# Path to your phenocodes file (should contain 'filename' and 'description' columns)
phenocodes_path <- "phenocodes" 

# Output heatmap filename
output_heatmap_filename <- "PGSscore_heatmap_Individuals_Normalized_AnnotatedRaw_SU_CroppedPheno.png" 

# Trait selection parameters (as in Python script)
PERCENTILE_THRESHOLD <- 99 # Filter by top N percentile of overall absolute differences across all individuals/traits
MAX_TRAITS_TO_PLOT <- 30 # Keep a cap for readability

# --- HELPER FUNCTION TO EXTRACT INDIVIDUAL ID FROM FILENAME (Python-style) ---
# Extracts e.g., "MC_M8_C-O" from "PGSscores_difference_MC_M8_C-O.txt"
extract_individual_id_python <- function(filepath) {
  filename <- basename(filepath)
  id <- str_replace(filename, "^PGSscores_difference_", "")
  id <- str_replace(id, "\\.txt$", "")
  return(id)
}

# --- HELPER FUNCTION TO CLEAN TRAIT LABELS (Python-style) ---
# Cleans e.g., "UKB_B_N3_FSH-both_sexes_filtered.tsv.bgz" to "UKB_B_N3_FSH"
clean_label_python <- function(label) {
  label_str <- as.character(label)
  label_str <- str_replace(label_str, "_filtered\\.tsv\\.bgz", "")
  # The split function needs careful handling in R compared to Python's [0]
  parts <- str_split(label_str, "-both_sexes", simplify = TRUE)
  return(parts[, 1])
}

# --- NEW HELPER FUNCTION: Crop string to the Nth underscore ---
crop_to_nth_underscore <- function(text, n) {
  if (is.na(text)) return(NA) # Handle NA inputs gracefully
  parts <- str_split(text, "_")[[1]]
  if (length(parts) <= n) {
    return(text) # If less than or equal to n underscores, return original string
  } else {
    return(paste(parts[1:n], collapse = "_")) # Recombine parts up to the nth
  }
}

# --- MAIN DATA PROCESSING FUNCTION FOR AN INDIVIDUAL FILE ---
# Reads an individual's score difference file, cleans labels, and merges phenocodes.
process_individual_file_python <- function(filepath, phenocodes_df, verbose = TRUE) {
  
  individual_id <- extract_individual_id_python(filepath)
  
  tryCatch({
    df_individual <- read_tsv(filepath, show_col_types = FALSE)
    
    required_cols <- c('match_file', 'score_difference')
    if (!all(required_cols %in% colnames(df_individual))) {
      stop(paste0("Missing required columns (", paste(required_cols, collapse = ", "), ") in ", filepath,
                  "\n  Found columns: ", paste(colnames(df_individual), collapse = ", ")))
    }
    
    df_processed <- df_individual %>%
      select(match_file = !!sym("match_file"), score_difference = !!sym("score_difference")) %>%
      filter(!is.na(score_difference) & !is.infinite(score_difference) & !is.nan(score_difference))
    
    if (nrow(df_processed) == 0) {
      if (verbose) message(paste("  - Warning: No valid score_difference data in '", filepath, "' after cleaning. Skipping."))
      return(NULL)
    }
    
    # Generate the base cleaned label (e.g., "UKB_B_N3_FSH")
    df_processed <- df_processed %>%
      mutate(
        label_clean_base = clean_label_python(match_file), # The code part (e.g., UKB_B_N3_FSH_v2)
        # Apply the new cropping rule to the base label for the code part of the final Phenotype_Label
        cropped_code_label = sapply(label_clean_base, crop_to_nth_underscore, n = 3),
        filename_clean = str_replace(match_file, "_filtered\\.tsv\\.bgz", ".tsv.bgz")
      )
    
    # Merge with phenocodes for descriptive labels
    # CORRECTED LINE: Specify that filename_clean from df_processed maps to filename in phenocodes_df
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
    
    # Handle duplicate trait labels by taking mean score_difference (as in Python)
    if (any(duplicated(merged_df$Phenotype_Label))) {
      if (verbose) message(paste0("  - Warning: Duplicate trait labels found in ", filepath, ". Taking mean score_difference for duplicates."))
      merged_df <- merged_df %>%
        group_by(Phenotype_Label) %>%
        summarise(score_difference = mean(score_difference, na.rm = TRUE), .groups = 'drop')
    }
    
    # Rename score_difference column to the individual_id for combining (similar to Python's .rename)
    df_final <- merged_df %>%
      select(Phenotype_Label, !!individual_id := score_difference) # Dynamic column naming
    
    return(df_final)
    
  }, error = function(e) {
    if (verbose) message(paste("  - An error occurred while reading individual file", filepath, ":", e$message))
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

# --- PROCESS ALL INDIVIDUAL FILES ---
all_individuals_dfs <- list()
message(paste0("\n--- Processing specified individual files ---"))

if (length(input_files) == 0) {
  stop("Error: No input files specified in the 'input_files' list. Exiting.")
}

for (filepath in input_files) {
  if (file.exists(filepath)) {
    message(paste0("  - Processing '", filepath, "'..."))
    processed_df <- process_individual_file_python(filepath, phenos)
    if (!is.null(processed_df)) {
      all_individuals_dfs[[length(all_individuals_dfs) + 1]] <- processed_df
    }
  } else {
    message(paste0("Warning: Input file not found: '", filepath, "'. Skipping."))
  }
}

if (length(all_individuals_dfs) == 0) {
  stop("No valid individual data found from input files. Exiting.")
}

# Combine all individual DataFrames into a single heatmap-ready DataFrame
# Traits as rows, individual IDs as columns (using full_join to align by Phenotype_Label)
message("\n--- Combining raw data from all individuals for heatmap ---")
heatmap_data_raw <- all_individuals_dfs %>%
  purrr::reduce(full_join, by = "Phenotype_Label") %>%
  rename(Phenotype = Phenotype_Label) # Rename for plotting consistency

message(paste0("  - Final heatmap_data_raw DataFrame shape (before filtering): ", nrow(heatmap_data_raw), "x", ncol(heatmap_data_raw)))


# --- TRAIT SELECTION (Percentile Filtering based on overall max absolute difference) ---
message("\n--- Trait Selection (Percentile Filtering) ---")

# Identify numeric columns (individual score difference columns)
numeric_cols <- names(heatmap_data_raw)[sapply(heatmap_data_raw, is.numeric)]

if (length(numeric_cols) == 0) {
  stop("No numeric score columns found after combining data. Check input files and processing.")
}

# Calculate the maximum absolute score difference for each trait across all individuals
max_abs_diff_per_trait_df <- heatmap_data_raw %>%
  select(Phenotype, all_of(numeric_cols)) %>%
  mutate(across(all_of(numeric_cols), ~replace(., is.infinite(.), NA))) %>%
  rowwise() %>% 
  mutate(max_abs_diff = max(abs(c_across(all_of(numeric_cols))), na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(desc(max_abs_diff)) %>%
  select(Phenotype, max_abs_diff)

# Remove any NaNs from the max_abs_diff_per_trait for quantile calculation
max_abs_diff_clean <- max_abs_diff_per_trait_df$max_abs_diff[!is.na(max_abs_diff_per_trait_df$max_abs_diff)]
names(max_abs_diff_clean) <- max_abs_diff_per_trait_df$Phenotype[!is.na(max_abs_diff_per_trait_df$max_abs_diff)]

total_unique_traits <- length(max_abs_diff_clean)

selected_labels <- character(0) # Initialize as empty character vector

if (total_unique_traits > 0) {
  percentile_value <- quantile(max_abs_diff_clean, 1 - (PERCENTILE_THRESHOLD / 100.0), na.rm = TRUE)
  
  if (percentile_value == 0.0) {
    qualifying_traits <- names(max_abs_diff_clean[max_abs_diff_clean > 0])
  } else {
    qualifying_traits <- names(max_abs_diff_clean[max_abs_diff_clean >= percentile_value])
  }
  
  # Ensure the selected labels respect MAX_TRAITS_TO_PLOT and are sorted by max_abs_diff
  selected_labels <- head(qualifying_traits[order(max_abs_diff_clean[qualifying_traits], decreasing = TRUE)], MAX_TRAITS_TO_PLOT)
  
  message(paste0("  - Total unique traits available: ", total_unique_traits))
  message(paste0("  - Threshold value for ", PERCENTILE_THRESHOLD, "th percentile (absolute difference): ", round(percentile_value, 4)))
  message(paste0("  - Number of traits meeting or exceeding the ", PERCENTILE_THRESHOLD, "th percentile: ", length(qualifying_traits)))
  
  if (length(selected_labels) == 0 && total_unique_traits > 0) {
    message(paste0("    (No traits strictly above ", round(percentile_value, 4), " or MAX_TRAITS_TO_PLOT is 0. Selecting top 1 trait as fallback.)"))
    selected_labels <- head(names(max_abs_diff_clean[order(max_abs_diff_clean, decreasing = TRUE)]), 1)
  }
} else {
  message("No unique traits available for selection criteria.")
}

message(paste0("  - Final number of labels selected for plotting (max ", MAX_TRAITS_TO_PLOT, "): ", length(selected_labels)))

if (length(selected_labels) == 0) {
  stop("No traits met the selection criteria for plotting. Cannot create heatmap. Exiting.")
}

# Filter the raw data to include only these selected labels for both color and annotation
filtered_heatmap_data_raw_for_annot <- heatmap_data_raw %>%
  filter(Phenotype %in% selected_labels)

# --- Z-SCORE NORMALIZATION FOR PLOTTING (Within Each Individual - Column-wise) ---
message("\n--- Normalizing score differences within each individual for heatmap ---")

normalized_heatmap_data <- filtered_heatmap_data_raw_for_annot %>%
  select(-Phenotype) # Exclude Phenotype column for normalization

# Apply Z-score normalization column-wise using scale() for each individual
normalized_heatmap_data <- as.data.frame(lapply(normalized_heatmap_data, function(col) {
  if (sd(col, na.rm = TRUE) == 0) {
    return(rep(0, length(col))) # Return a vector of zeros if SD is zero
  } else {
    return(scale(col, center = TRUE, scale = TRUE)) # Standard Z-score scaling
  }
}))
# Reassign column names lost by as.data.frame, exclude Phenotype column that was dropped
colnames(normalized_heatmap_data) <- names(filtered_heatmap_data_raw_for_annot)[-1] 
normalized_heatmap_data$Phenotype <- filtered_heatmap_data_raw_for_annot$Phenotype # Add Phenotype column back

# Convert to long format for ggplot for normalized data
normalized_heatmap_data_long <- normalized_heatmap_data %>%
  pivot_longer(
    cols = -Phenotype,
    names_to = "SampleID",
    values_to = "Z_Score"
  )

# Prepare annotation data (raw values) in long format and ensure alignment
filtered_heatmap_data_raw_for_annot_long <- filtered_heatmap_data_raw_for_annot %>%
  pivot_longer(
    cols = -Phenotype,
    names_to = "SampleID",
    values_to = "Raw_Score_Difference"
  )

# Join them so we have Z_Score (for color) and Raw_Score_Difference (for text) in one dataframe
plot_data <- normalized_heatmap_data_long %>%
  left_join(filtered_heatmap_data_raw_for_annot_long, by = c("Phenotype", "SampleID"))

if (nrow(plot_data) == 0) {
  stop("No data remaining after normalization and joining. Cannot create heatmap. Exiting.")
}

# Sort traits by their overall normalized score for consistent plotting order
# Calculate mean of Z_Score for each Phenotype, then order the factor levels
phenotype_mean_z <- plot_data %>%
  group_by(Phenotype) %>%
  summarise(overall_normalized_score = mean(Z_Score, na.rm = TRUE)) %>%
  arrange(desc(overall_normalized_score))

plot_data$Phenotype <- factor(plot_data$Phenotype, levels = phenotype_mean_z$Phenotype)

# Order SampleID factor levels based on their original appearance in input_files
sample_order_from_files <- sapply(input_files, extract_individual_id_python)
plot_data$SampleID <- factor(plot_data$SampleID, levels = sample_order_from_files)


message(paste0("  - Filtered and Normalized heatmap_data DataFrame (long format) shape: ", nrow(plot_data), "x", ncol(plot_data)))


# --- PLOTTING THE HEATMAP ---
message("\n--- Generating heatmap with normalized individual differences ---")

# Dynamic sizing based on number of traits and individuals for readability
num_phenotypes <- n_distinct(plot_data$Phenotype)
num_samples <- n_distinct(plot_data$SampleID)

base_width <- 8 # Minimum width
base_height <- 6 # Minimum height

# Scale dimensions based on number of rows/columns
plot_width_in <- max(base_width, num_samples * 1.5) # Adjust cell width for readability
plot_height_in <- max(base_height, num_phenotypes * 0.3) # Adjust cell height for readability

# Cap max height/width to avoid excessively large plots (can be adjusted)
if (plot_height_in > 30) plot_height_in <- 30 
if (plot_width_in > 20) plot_width_in <- 20


p_heatmap <- ggplot(plot_data, aes(x = SampleID, y = Phenotype, fill = Z_Score)) +
  geom_tile(color = "black", linewidth = 0.5) + # Tiles with black borders (like Python's linecolor)
  
  # Annotate with raw score differences, formatted to 0 decimal places
  #geom_text(aes(label = round(Raw_Score_Difference, 0)), 
   #         color = "black", # Text color for annotations
    #        size = 15 / .pt, # Equivalent to fontsize=15 in matplotlib (convert mm to pts)
     #       na.rm = TRUE) + # Don't plot text for NA values
  
  scale_fill_gradient2(
    low = "blue",      # Blue for negative (decrease)
    mid = "white",     # White for mid (near zero)
    high = "red",      # Red for positive (increase)
    midpoint = 0,      # Center the color scale at 0
    name = "Normalized Score Difference\n(Z-score per individual)" # Colorbar label
  ) +
  labs(
    title = paste0("Heatmap of PGS Score Differences per Individual\n(Colors: Normalized Z-score- Top ", 
                   num_phenotypes, " Traits by ", PERCENTILE_THRESHOLD, "th Percentile Absolute Change)"),
    x = "Individual ID",
    y = "PGS Trait"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray30"), 
    axis.title.x = element_text(size = 18, face = "bold", color = "black"), # Larger font
    axis.title.y = element_text(size = 18, face = "bold", color = "black"), # Larger font
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 16, color = "black", face = "bold"), # Rotate x-axis labels, larger font, bold
    axis.text.y = element_text(size = 10, color = "black", face = "bold"), # Bold, smaller font for y-axis
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    panel.grid.major = element_blank(), # Remove grid lines within the heatmap
    panel.grid.minor = element_blank(),
    axis.ticks = element_blank()
  )

# --- SAVE THE PLOT ---
ggsave(output_heatmap_filename,
       plot = p_heatmap,
       width = plot_width_in,
       height = plot_height_in,
       units = "in", dpi = 300, bg = "white") # dpi 300 as in Python

message(paste0("\n--- Heatmap saved to '", output_heatmap_filename, "' ---"))

# To display the plot
print(p_heatmap)

message("--- Plotting complete ---")