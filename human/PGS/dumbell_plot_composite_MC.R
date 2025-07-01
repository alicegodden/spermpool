# remove behavioural terms
# PGS Score Change Scatter Plots (Dumbbell Plots) for Paired Raw Scores - Composite Plot

# --- INSTALL & LOAD NECESSARY LIBRARIES ---
# Ensure patchwork is installed for composite plotting
if (!requireNamespace("patchwork", quietly = TRUE)) {
  install.packages("patchwork")
}
library(dplyr)      # For data manipulation
library(readr)      # For reading data (read_tsv)
library(ggplot2)    # For plotting
library(stringr)    # For string manipulation
library(tools)      # For file path manipulation (basename)
library(tidyr)      # For pivot_longer, pivot_wider
library(purrr)      # For reduce (to combine data frames efficiently)
library(patchwork)  # For combining plots

# --- SET WORKING DIRECTORY ---
# <<< IMPORTANT: Adjust this path to where your input files and phenocodes file are located >>>
setwd("~/PycharmProjects/polygenic/PGS_scores_test/SU/T0") 

# --- CONFIGURATION ---
# >>> IMPORTANT: Adjust these to match your actual file structure and desired behavior <<<

# Input files (raw score .tsv.gz files, assumed to be paired C, O, C, O...)
input_files <- c(
  "GRCh37_M8_C1.dedup_bcftools__CHR_nc_to_chr.nochr.aaf_scores.tsv.gz",
  "GRCh37_M8_O1.dedup_bcftools__CHR_nc_to_chr.nochr.aaf_scores.tsv.gz",
  "GRCh37_M11_C2._merged_dedup_bcftools__CHR_nc_to_chr.nochr.aaf_scores.tsv.gz", # Corrected M11 C2 filename based on typical pattern
  "GRCh37_M11_O1.dedup_bcftools__CHR_nc_to_chr.nochr.aaf_scores.tsv.gz",
  "GRCh37_M12_C1.dedup_bcftools__CHR_nc_to_chr.nochr.aaf_scores.tsv.gz",
  "GRCh37_M12_O1.dedup_bcftools__CHR_nc_to_chr.nochr.aaf_scores.tsv.gz",
  "GRCh37_M13_C1.dedup_bcftools__CHR_nc_to_chr.nochr.aaf_scores.tsv.gz",
  "GRCh37_M13_O1.dedup_bcftools__CHR_nc_to_chr.nochr.aaf_scores.tsv.gz"
)

# Path to your phenocodes file (should contain 'filename' and 'description' columns)
phenocodes_path <- "phenocodes" 

# Output composite plot filename
output_composite_plot_filename <- "PGS_Change_Dumbbell_Composite_PerPairTopPhenotypes.png" 

# Trait selection parameters
PERCENTILE_THRESHOLD <- 99 # Filter by top N percentile of absolute differences (now per pair)
MAX_TRAITS_TO_PLOT <- 15   # Cap for readability

# --- CONFIGURABLE COLUMN NAMES IN YOUR INPUT TSV.GZ FILES ---
# !!! IMPORTANT: Adjust these to match the actual column names in your files !!!
PHENOTYPE_COL_NAME <- "match_file" 
SCORE_COL_NAME <- "product_total"  

# --- HELPER FUNCTION TO EXTRACT INDIVIDUAL ID FROM FILENAME ---
# Extracts the segment between the first underscore and the first dot.
extract_individual_id <- function(filepath) {
  filename <- basename(filepath)
  # 1. Get the part of the string after the first underscore
  after_first_underscore <- str_replace(filename, "^[^_]*_", "")
  # 2. From this new string, get the part before the first dot
  id <- str_replace(after_first_underscore, "\\..*$", "")
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
# 'fixed = TRUE' is removed here because we are using a regex pattern with '|'
terms_to_remove_pattern <- "Place of birth|television|Nap|Getting up|Z34|Alcohol|Time spent using computer|Hot drink|Lifetime number of sexual partners|fruit|Lamb|Worry|Cereal|PM: final answer"
all_raw_data_long <- all_raw_data_long %>%
  filter(!str_detect(Phenotype_Label, terms_to_remove_pattern))

filtered_pheno_count <- length(unique(all_raw_data_long$Phenotype_Label))
message(paste0("  - Removed phenotypes containing: '", gsub("\\|", "', '", terms_to_remove_pattern), "'. Initial unique phenotypes: ", initial_pheno_count, ", After filter: ", filtered_pheno_count, "."))


#message(paste0("  - Total raw data rows loaded: ", nrow(all_raw_data_long)))
# --- Generate Plots for Each Pair and Store in a List ---
message("\n--- Generating individual dumbbell plots for each donor pair ---")
plot_list <- list()

for (i in seq(1, length(input_files), by = 2)) {
  file_c_path <- input_files[i]
  file_o_path <- input_files[i+1]
  
  sample_c_id <- extract_individual_id(file_c_path)
  sample_o_id <- extract_individual_id(file_o_path)
  donor_id <- str_extract(sample_c_id, "^M\\d+") # Extract donor ID (e.g., M8)
  
  message(paste0("  - Preparing plot for donor pair: ", donor_id, " (", sample_c_id, " vs ", sample_o_id, ")"))
  
  # Filter data for the current pair
  current_pair_data <- all_raw_data_long %>%
    filter(SampleID %in% c(sample_c_id, sample_o_id)) %>%
    pivot_wider(names_from = SampleID, values_from = Score) %>%
    # Filter out rows where either C or O score is missing for the selected phenotypes
    filter(!is.na(.data[[sample_c_id]]) & !is.na(.data[[sample_o_id]])) %>%
    mutate(
      Score_C = .data[[sample_c_id]],
      Score_O = .data[[sample_o_id]],
      absolute_difference = abs(Score_O - Score_C)
    )
  
  if (nrow(current_pair_data) == 0) {
    message(paste0("    - No complete data for donor pair ", donor_id, ". Skipping plot."))
    # Add an empty plot or a placeholder for missing data to maintain layout
    plot_list[[paste0("plot_", donor_id)]] <- ggplot() + 
      labs(title = paste0(donor_id, " (No Data)"), 
           x = "Raw PGS Score", y = "Phenotype") + 
      theme_void() + 
      theme(plot.title = element_text(hjust = 0.5, size = 12, face = "italic"))
    next
  }
  
  # --- Per-Pair Trait Selection (NEW LOGIC) ---
  message(paste0("    - Selecting top phenotypes for pair ", donor_id, " (", PERCENTILE_THRESHOLD, "th percentile)..."))
  
  # Determine the 99th percentile threshold for absolute differences in this pair
  percentile_value_per_pair <- quantile(current_pair_data$absolute_difference, 1 - (PERCENTILE_THRESHOLD / 100.0), na.rm = TRUE)
  
  # Select top phenotypes for this specific pair
  selected_top_phenotypes_for_this_pair <- current_pair_data %>%
    filter(absolute_difference >= percentile_value_per_pair) %>%
    arrange(desc(absolute_difference)) %>% # Order by largest difference
    head(MAX_TRAITS_TO_PLOT) %>%
    pull(Phenotype_Label)
  
  if (length(selected_top_phenotypes_for_this_pair) == 0) {
    message(paste0("    - No phenotypes met the ", PERCENTILE_THRESHOLD, "th percentile threshold for pair ", donor_id, ". Selecting top 1 trait as fallback."))
    selected_top_phenotypes_for_this_pair <- head(current_pair_data %>% arrange(desc(absolute_difference)) %>% pull(Phenotype_Label), 1)
  }
  
  message(paste0("    - Selected ", length(selected_top_phenotypes_for_this_pair), " top phenotypes for pair ", donor_id, "."))
  
  # Filter data for the current plot based on the *per-pair* selected phenotypes
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
  
  # Order phenotypes on Y-axis by their 'Control' score for this pair
  phenotype_order_for_plot <- pair_plot_data_filtered %>%
    arrange(Score_C) %>%
    pull(Phenotype_Label)
  
  pair_plot_data_filtered$Phenotype_Label <- factor(pair_plot_data_filtered$Phenotype_Label, levels = phenotype_order_for_plot)
  
  # Prepare data for plotting points (long format)
  plot_points_data <- pair_plot_data_filtered %>%
    pivot_longer(
      cols = c(Score_C, Score_O),
      names_to = "Score_Type",
      values_to = "Score_Value"
    ) %>%
    mutate(
      Dot_Color_Category = case_when( # This will drive the color aesthetic
        Score_Type == "Score_C" ~ "Control",
        Score_Type == "Score_O" ~ "Outcome"
      )
    )
  
  # Determine X-axis limits for this specific plot
  plot_min_score <- min(plot_points_data$Score_Value, na.rm = TRUE)
  plot_max_score <- max(plot_points_data$Score_Value, na.rm = TRUE)
  
  # Add padding to x-axis
  padding <- (plot_max_score - plot_min_score) * 0.1
  if (padding == 0) padding <- 1 # Fallback for very narrow ranges
  plot_min_score_padded <- plot_min_score - padding
  plot_max_score_padded <- plot_max_score + padding
  
  # Ensure min/max are not identical after padding
  if (plot_min_score_padded == plot_max_score_padded) {
    plot_min_score_padded <- plot_min_score - 0.1
    plot_max_score_padded <- plot_max_score + 0.1
  }
  
  p <- ggplot(pair_plot_data_filtered, aes(y = Phenotype_Label)) +
    # Draw lines connecting the points (fixed neutral color)
    geom_segment(aes(x = Score_C, xend = Score_O), 
                 linewidth = 0.8, color = "grey50") +
    
    # Draw points for both Control and Outcome, colored by type, AND SHAPED BY TYPE
    geom_point(data = plot_points_data, 
               aes(x = Score_Value, color = Dot_Color_Category, shape = Dot_Color_Category), # <-- MODIFIED LINE: Added shape = Dot_Color_Category
               size = 3, alpha=0.6) + # <-- Removed fixed shape = 16
    
    labs(
      subtitle = paste0("Donor: ", donor_id), # Subtitle for individual plot
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
    # NEW: Add scale_shape_manual to define shapes for Control and Outcome
    scale_shape_manual(
      name = "Sample Type", # Use the same legend title as scale_color_manual for a combined legend
      values = c("Control" = 16, "Outcome" = 17), # 16 is a filled circle, 17 is a filled triangle
      labels = c(paste0("Control (", sample_c_id, ")"), paste0("Outcome (", sample_o_id, ")")) # Match labels for consistency
    )
  
  # Store plot in list
  plot_list[[paste0("plot_", donor_id)]] <- p
}

# --- Create Composite Plot ---
message("\n--- Creating composite plot ---")

# Determine optimal layout for patchwork (e.g., 2 columns)
n_plots <- length(plot_list)
n_cols_composite <- 2 # Number of columns for the composite plot
n_rows_composite <- ceiling(n_plots / n_cols_composite)

# Combine plots using patchwork
# Use `wrap_plots` to arrange the plots
composite_plot <- wrap_plots(plot_list, ncol = n_cols_composite) +
  plot_annotation(
    title = paste0("PGS Raw Score Changes by donor\n(Top ", MAX_TRAITS_TO_PLOT, " Phenotypes per Pair, displaying 99th percentile)"),
    theme = theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold", margin = margin(b = 10))) # Overall title style
  ) & theme(legend.position = "bottom") # Apply legend position to all subplots

# Calculate dynamic height for the composite plot
# Base height per row, plus some extra per trait for long plots
height_per_panel <- max(7, 0.25 * MAX_TRAITS_TO_PLOT) # MAX_TRAITS_TO_PLOT is now fixed for height
composite_plot_height_in <- n_rows_composite * height_per_panel
composite_plot_width_in <- n_cols_composite * 8 # Assuming 8 inches base width per panel

# Cap max height/width to avoid excessively large plots
if (composite_plot_height_in > 40) composite_plot_height_in <- 40 
if (composite_plot_width_in > 30) composite_plot_width_in <- 30 

# --- SAVE THE COMPOSITE PLOT ---
ggsave(output_composite_plot_filename,
       plot = composite_plot,
       width = composite_plot_width_in,
       height = composite_plot_height_in,
       units = "in", dpi = 300, bg = "white")

message(paste0("\n--- Composite plot saved to '", output_composite_plot_filename, "' ---"))

# To display the plot (might be very large if many panels/traits)
print(composite_plot)

message("--- Plotting complete ---")