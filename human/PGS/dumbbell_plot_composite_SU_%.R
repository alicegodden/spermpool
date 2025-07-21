Title: dumbbell_plot_composite_SU_%.R
Author: Dr. Alice M. Godden

## custom x axes for standardisation

# --- SCRIPT METADATA ---
# Project: PGS Score Change Visualization
# Description: Generates dumbbell plots to visualize PGS raw score changes for paired
#              samples (Control vs. Outcome). Features include fixed X-axis range (0-0.5),
#              Y-axis reordering, and explicit display of percentage change.
#              Phenotype labels are now directly on the plot.
# Author: [Your Name/Org - if applicable]
# Date: July 17, 2025 (Updated) - Corrected as of July 21, 2025

# --- INSTALL & LOAD NECESSARY LIBRARIES ---
# Ensure patchwork is installed for composite plotting
if (!requireNamespace("patchwork", quietly = TRUE)) {
  install.packages("patchwork")
}
# Optional: Install 'here' for portable path management
if (!requireNamespace("here", quietly = TRUE)) {
  install.packages("here")
}

library(dplyr)       # For data manipulation (e.g., select, filter, mutate, group_by, summarise)
library(readr)       # For reading data (read_tsv)
library(ggplot2)     # For plotting (geom_segment, geom_point, theme_minimal, labs, etc.)
library(stringr)     # For string manipulation (str_replace, str_split, str_detect, str_extract)
library(tools)       # For file path manipulation (basename)
library(tidyr)       # For data reshaping (pivot_longer, pivot_wider)
library(purrr)       # For functional programming (reduce - though not explicitly used in final iteration, good to keep)
library(patchwork)   # For combining multiple ggplot objects into a composite plot
library(here)        # For relative file paths (optional, uncomment usage below)


# --- SET WORKING DIRECTORY / PATHS ---
# <<< IMPORTANT: Adjust this path to where your input files and phenocodes file are located >>>
# Using 'here' package for portable path management. Uncomment and adjust if preferred.
# setwd(here::here("PycharmProjects", "polygenic", "PGS_scores_test", "SU", "T0"))
# If not using 'here', use direct path:
setwd("~/PycharmProjects/polygenic/PGS_scores_test/SU/T0")


# --- CONFIGURATION PARAMETERS ---
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

# Output composite plot filename
output_composite_plot_filename <- "PGS_PercentageChange_Dumbbell_FixedAxis_Composite.png"

# Trait selection parameters
PERCENTILE_THRESHOLD <- 99 # Filter by top N percentile of absolute differences (now per pair)
MAX_TRAITS_TO_PLOT <- 15    # Cap for readability (maximum phenotypes shown per donor)

# --- FIXED X-Axis Standardization Parameters (Applied Globally) ---
X_AXIS_MIN <- -0.5
X_AXIS_MAX <- 0.5

# --- CONFIGURABLE COLUMN NAMES IN YOUR INPUT TSV.GZ FILES ---
# !!! IMPORTANT: Adjust these to match the actual column names in your files !!!
PHENOTYPE_COL_NAME <- "match_file"
SCORE_COL_NAME <- "product_total"

# --- HELPER FUNCTIONS ---

#' Extracts the individual ID from a given file path.
#' Updated to extract D#T# pattern.
extract_individual_id <- function(filepath) {
  filename <- basename(filepath)
  # Extract "D#T#" or "D#T#a/b" pattern
  id <- str_extract(filename, "D\\d+T\\d+[ab]?") # Uses [ab]? to capture 'a' or 'b'
  return(id)
}

#' Cleans trait labels, removing specific suffixes and 'both_sexes' markers.
clean_label_python <- function(label) {
  label_str <- as.character(label)
  label_str <- str_replace(label_str, "_filtered\\.tsv\\.bgz", "")
  parts <- str_split(label_str, "-both_sexes", simplify = TRUE)
  return(parts[, 1])
}

#' Crops a string to the Nth underscore.
#' Useful for creating shorter, more readable phenotype codes.
crop_to_nth_underscore <- function(text, n) {
  if (is.na(text)) return(NA)
  parts <- str_split(text, "_")[[1]]
  if (length(parts) <= n) {
    return(text)
  } else {
    return(paste(parts[1:n], collapse = "_"))
  }
}

#' Creates a placeholder ggplot object.
create_placeholder_plot <- function(donor_id, message_text = "No Data") {
  ggplot(data.frame(x=0, y=0)) +
    geom_point(aes(x,y), color="transparent") + # Invisible point to keep scale
    labs(title = paste0(donor_id, " (", message_text, ")"),
         x = "PGS Score Change", y = "Phenotype") +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 12, face = "italic"),
      axis.title.x = element_text(size = 10, face = "bold", color = "black"),
      axis.title.y = element_text(size = 10, face = "bold", color = "black"),
      plot.margin = unit(c(0.1, 0.25, 0.1, 0.25), "cm")
    )
}


#' Processes a single raw score file.
#' Reads the .tsv.gz, cleans phenotype labels, handles duplicates, and attaches SampleID.
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
        # Construct the final Phenotype_Label: only description if available, otherwise cropped code
        Phenotype_Label = ifelse(
          !is.na(description),
          description, # Only the description
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

# --- MAIN SCRIPT EXECUTION ---

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

message(paste0("  - Total raw data rows loaded: ", nrow(all_raw_data_long)))

# --- Generate Plots for Each Pair and Store in a List ---
message("\n--- Generating individual dumbbell plots for each donor pair ---")
plot_list <- list()

for (i in seq(1, length(input_files), by = 2)) {
  file_c_path <- input_files[i]
  file_o_path <- input_files[i + 1]
  
  sample_c_id <- extract_individual_id(file_c_path)
  sample_o_id <- extract_individual_id(file_o_path)
  
  # --- FIX: More robust donor_id extraction for plot keys and custom chunks ---
  # Extract the "D#" part
  d_part <- str_extract(sample_c_id, "D\\d+")
  # Extract the optional 'a' or 'b' part specifically
  ab_part <- str_extract(sample_c_id, "[ab]") # Use [ab] to explicitly look for 'a' or 'b'
  
  # Combine them to form the final donor_id (e.g., D6a, D6b, D1, D2)
  donor_id <- if (!is.na(ab_part)) {
    paste0(d_part, ab_part)
  } else {
    d_part
  }
  # --- END FIX ---
  
  message(paste0("  - DEBUG: i=", i, ", file_c_path='", file_c_path, "', sample_c_id='", sample_c_id, "', derived donor_id='", donor_id, "'"))
  
  message(paste0("  - Preparing plot for donor pair: ", donor_id, " (", sample_c_id, " vs ", sample_o_id, ")"))
  
  # Filter data for the current pair BEFORE pivoting
  current_pair_raw_filtered <- all_raw_data_long %>%
    filter(SampleID %in% c(sample_c_id, sample_o_id))
  
  # --- CRITICAL FIX: Add robustness check here for missing sample data ---
  # Check if both sample IDs are actually present in the filtered data.
  if (length(unique(current_pair_raw_filtered$SampleID)) < 2) {
    message(paste0("    - Warning: Incomplete data for donor pair ", donor_id, ". Both samples (", sample_c_id, ", ", sample_o_id, ") not found. Adding placeholder."))
    plot_list[[paste0("plot_", donor_id)]] <- create_placeholder_plot(donor_id, "Incomplete Pair Data")
    message(paste0("  - DEBUG: Plot key '", paste0("plot_", donor_id), "' added/updated."))
    next
  }
  # --- END CRITICAL FIX ---
  
  current_pair_data <- current_pair_raw_filtered %>%
    pivot_wider(names_from = SampleID, values_from = Score) %>%
    # Filter out rows where either C or O score is missing for the selected phenotypes
    filter(!is.na(.data[[sample_c_id]]) & !is.na(.data[[sample_o_id]])) %>%
    mutate(
      Score_C = .data[[sample_c_id]],
      Score_O = .data[[sample_o_id]],
      absolute_difference = abs(Score_O - Score_C) # Keep for filtering
    )
  
  if (nrow(current_pair_data) == 0) {
    message(paste0("    - No complete data for donor pair ", donor_id, " after NA filtering. Adding placeholder."))
    plot_list[[paste0("plot_", donor_id)]] <- create_placeholder_plot(donor_id, "No Valid Pairs")
    message(paste0("  - DEBUG: Plot key '", paste0("plot_", donor_id), "' added/updated."))
    next
  }
  
  # --- Per-Pair Trait Selection ---
  message(paste0("    - Selecting top phenotypes for pair ", donor_id, " (", PERCENTILE_THRESHOLD, "th percentile)..."))
  
  percentile_value_per_pair <- quantile(current_pair_data$absolute_difference, 1 - (PERCENTILE_THRESHOLD / 100.0), na.rm = TRUE)
  
  selected_top_phenotypes_df <- current_pair_data %>%
    filter(absolute_difference >= percentile_value_per_pair) %>%
    arrange(desc(absolute_difference)) %>% # Order by largest difference
    head(MAX_TRAITS_TO_PLOT) # Cap the number of traits to plot
  
  if (nrow(selected_top_phenotypes_df) == 0) {
    message(paste0("    - No phenotypes met the ", PERCENTILE_THRESHOLD, "th percentile threshold for pair ", donor_id, ". Selecting top 1 trait as fallback, if available."))
    # Fallback to top 1 if no percentile data, otherwise an empty placeholder
    if(nrow(current_pair_data) > 0) { # Check if there was any data before percentile filter
      selected_top_phenotypes_df <- head(current_pair_data %>% arrange(desc(absolute_difference)), 1)
      message(paste0("      - Fallback: Plotting top 1 phenotype for ", donor_id, "."))
    } else {
      message(paste0("      - No data available even for fallback for ", donor_id, ". Adding placeholder."))
      plot_list[[paste0("plot_", donor_id)]] <- create_placeholder_plot(donor_id, "No Top Phenos")
      message(paste0("  - DEBUG: Plot key '", paste0("plot_", donor_id), "' added/updated."))
      next
    }
  }
  
  # --- Calculate Percentage Change ---
  selected_top_phenotypes_df <- selected_top_phenotypes_df %>%
    mutate(
      Score_C_perc_change = ifelse(Score_C == 0, 0, (Score_C - Score_C) / Score_C * 100), # Should always be 0
      Score_O_perc_change = ifelse(Score_C == 0, 0, (Score_O - Score_C) / Score_C * 100)
    )
  # Handle potential Inf/-Inf from division by zero if Score_C is 0 for some rows
  selected_top_phenotypes_df <- selected_top_phenotypes_df %>%
    mutate(
      Score_C_perc_change = replace(Score_C_perc_change, is.infinite(Score_C_perc_change) | is.nan(Score_C_perc_change), 0),
      Score_O_perc_change = replace(Score_O_perc_change, is.infinite(Score_O_perc_change) | is.nan(Score_O_perc_change), 0)
    )
  
  message(paste0("    - Selected ", nrow(selected_top_phenotypes_df), " top phenotypes for pair ", donor_id, "."))
  
  # Order phenotypes by their percentage change for consistent Y-axis ordering
  selected_top_phenotypes_df <- selected_top_phenotypes_df %>%
    arrange(Score_O_perc_change) %>% # Order by the percentage change of the outcome score
    mutate(Phenotype_Label = factor(Phenotype_Label, levels = unique(Phenotype_Label))) # Preserve order as factor levels
  
  # Prepare data for plotting points (long format)
  plot_points_data <- selected_top_phenotypes_df %>%
    pivot_longer(
      cols = c(Score_C_perc_change, Score_O_perc_change), # Use percentage change columns
      names_to = "Score_Type",
      values_to = "Score_Value"
    ) %>%
    mutate(
      Dot_Color_Category = case_when(
        Score_Type == "Score_C_perc_change" ~ "T0",
        Score_Type == "Score_O_perc_change" ~ "T4"
      ),
      # Add percentage change label directly for geom_text
      display_label = paste0(format(round(Score_Value, 1), nsmall = 1), "%")
    )
  
  # --- Create the plot for the current donor pair ---
  p_donor <- ggplot(selected_top_phenotypes_df, aes(y = Phenotype_Label)) +
    geom_segment(aes(x = Score_C_perc_change, xend = Score_O_perc_change), linewidth = 0.8, color = "grey50") +
    geom_point(data = plot_points_data,
               aes(x = Score_Value, color = Dot_Color_Category, shape = Dot_Color_Category),
               size = 3, alpha=0.9) +
    geom_text(data = plot_points_data %>% filter(Dot_Color_Category == "T4"), # Only label T4 points
              aes(x = Score_Value, label = display_label, color = Dot_Color_Category),
              hjust = -0.15, vjust = 0.5, size = 0, fontface = "bold", show.legend = FALSE) + # Adjust hjust to place label next to dot
    labs(
      title = paste0("Donor: ", donor_id), # Plot title for each donor
      x = "Percentage Change from T0 (%)", # Consistent X-axis label
      y = "Phenotype" # Explicitly set Y-axis title here
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
      axis.title.x = element_text(size = 10, face = "bold", color = "black"),
      axis.title.y = element_text(size = 10, face = "bold", color = "black"), # Ensure Y-axis title is visible
      axis.text.x = element_text(size = 8, color = "black", face = "bold"),
      axis.ticks.x = element_line(),
      axis.text.y = element_text(size = 10, color = "black", face = "bold"), # Y-axis label size
      legend.position = "none", # Legends will be collected by patchwork globally
      legend.title = element_text(size = 9, face = "bold"),
      legend.text = element_text(size = 8, face="bold"),
      panel.grid.major.y = element_line(color = "lightgray", linetype = "dotted", linewidth = 0.5),
      panel.grid.major.x = element_line(color = "lightgray", linetype = "dotted", linewidth = 0.5),
      panel.grid.minor = element_blank(),
      plot.margin = unit(c(0.1, 0.25, 0.1, 0.25), "cm")
    ) +
    coord_cartesian(xlim = c(X_AXIS_MIN, X_AXIS_MAX)) + # Fixed X-axis range
    scale_color_manual(values = c("T0" = "red", "T4" = "blue")) +
    scale_shape_manual(values = c("T0" = 16, "T4" = 17))
  
  plot_list[[paste0("plot_", donor_id)]] <- p_donor
  message(paste0("  - DEBUG: Plot key '", paste0("plot_", donor_id), "' added/updated."))
}

# --- DEBUGGING: Check plot_list keys after loop ---
message("\n--- Debugging: Keys in plot_list after generation ---")
print(names(plot_list))
message(paste0("  - Total plots generated (by unique keys): ", length(plot_list)))
# --- END DEBUGGING ---

# --- Create Composite Plot ---
message("\n--- Creating composite plot ---")

# Ensure plot_list is not empty before proceeding
if (length(plot_list) == 0) {
  stop("Error: 'plot_list' is empty. No plots were successfully generated for any donor pair. Cannot create composite plot.")
}

n_plots <- length(plot_list)
n_cols_composite <- 2
n_rows_composite <- ceiling(n_plots / n_cols_composite)

composite_plot <- wrap_plots(plot_list, ncol = n_cols_composite) +
  plot_annotation(
    title = paste0("PGS Percentage Change by Donor (Relative to T0)\nTop ", MAX_TRAITS_TO_PLOT, " Phenotypes (", PERCENTILE_THRESHOLD, "th percentile)"),
    theme = theme(plot.title = element_text(hjust = 0.5, size = 11, face = "bold", margin = margin(b = 10)))
  ) & theme(legend.position = "bottom")

# Calculate dynamic height for the composite plot (simplified for fixed axis)
# Each individual plot will have MAX_TRAITS_TO_PLOT phenotypes.
# We need to account for plot title, axes, and legend.
height_per_individual_plot <- (MAX_TRAITS_TO_PLOT * 0.25) + 2.5 # 0.25 inches per phenotype + 2.5 for overhead (title, axes, padding)

composite_plot_height_in <- n_rows_composite * height_per_individual_plot
composite_plot_height_in <- max(composite_plot_height_in, 8) # Ensure a reasonable minimum height
composite_plot_height_in <- min(composite_plot_height_in, 60) # Cap maximum height to prevent excessively large plots

composite_plot_width_in <- n_cols_composite * 8 # Adjusted width for fixed axis, 8 inches per column
if (composite_plot_width_in > 40) composite_plot_width_in <- 40 # Cap width

message(paste0("  - Saving composite plot to '", output_composite_plot_filename, "' with width=", round(composite_plot_width_in, 2), "in, height=", round(composite_plot_height_in, 2), "in..."))

tryCatch({
  ggsave(output_composite_plot_filename,
         plot = composite_plot,
         width = composite_plot_width_in,
         height = composite_plot_height_in,
         units = "in", dpi = 300, bg = "white")
  message(paste0("\n--- Script finished. Composite plot saved to: ", output_composite_plot_filename, " ---"))
}, error = function(e) {
  message(paste0("\n--- ERROR: Failed to save composite plot. Reason: ", e$message, " ---"))
  message("--- Please check the plot object and output path. ---")
})

# Print the composite plot object to display it in RStudio (if running interactively)
# This line might still cause issues in some RStudio versions for complex patchwork objects,
# even if ggsave succeeds. Comment out if the GGSAVE works and interactive display is not critical.
composite_plot
