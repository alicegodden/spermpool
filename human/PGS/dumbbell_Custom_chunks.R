# --- SCRIPT METADATA ---
# Project: PGS Score Change Visualization
# Description: Generates dumbbell plots to visualize PGS raw score changes for paired
#              samples (Control vs. Outcome). Features include Y-axis segmentation,
#              X-axis adjusted per chunk, Y-axis reordering, and custom chunk sizes
#              per donor to improve readability for top phenotypes.
# Author: Dr. Alice M. Godden


# --- INSTALL & LOAD NECESSARY LIBRARIES ---
# Ensure patchwork is installed for composite plotting
if (!requireNamespace("patchwork", quietly = TRUE)) {
  install.packages("patchwork")
}
# Optional: Install 'here' for portable path management
if (!requireNamespace("here", quietly = TRUE)) {
  install.packages("here")
}

library(dplyr)      # For data manipulation (e.g., select, filter, mutate, group_by, summarise)
library(readr)      # For reading data (read_tsv)
library(ggplot2)    # For plotting (geom_segment, geom_point, theme_minimal, labs, etc.)
library(stringr)    # For string manipulation (str_replace, str_split, str_detect, str_extract)
library(tools)      # For file path manipulation (basename)
library(tidyr)      # For data reshaping (pivot_longer, pivot_wider)
library(purrr)      # For functional programming (reduce - though not explicitly used in final iteration, good to keep)
library(patchwork)  # For combining multiple ggplot objects into a composite plot
library(here)       # For relative file paths (optional, uncomment usage below)


# --- SET WORKING DIRECTORY / PATHS ---
# <<< IMPORTANT: Adjust this path to where your input files and phenocodes file are located >>>
# Option 1: Set absolute path (less portable)
setwd("~/PycharmProjects/polygenic/PGS_scores_test/SU/T0")

# Option 2 (Recommended for portability): Use 'here' package
# Create a .here file in your project root (e.g., in PycharmProjects/polygenic/)
# and then define paths relative to that root.
# For example, if your .here file is in 'polygenic', and your data is in 'polygenic/PGS_scores_test/SU/T0':
# base_data_path <- here("PGS_scores_test", "SU", "T0")
# input_files_full_paths <- file.path(base_data_path, input_files)
# phenocodes_path_full <- file.path(base_data_path, phenocodes_path)
# For simplicity, sticking to setwd() as in your original script, but `here` is better for projects.


# --- CONFIGURATION PARAMETERS ---
# >>> IMPORTANT: Adjust these to match your actual file structure and desired behavior <<<

# Input files (raw score .tsv.gz files, assumed to be paired C, O, C, O...)
input_files <- c(
  "GRCh37_M8_C1.dedup_bcftools__CHR_nc_to_chr.nochr.aaf_scores.tsv.gz",
  "GRCh37_M8_O1.dedup_bcftools__CHR_nc_to_chr.nochr.aaf_scores.tsv.gz",
  "GRCh37_M11_C2._merged_dedup_bcftools__CHR_nc_to_chr.nochr.aaf_scores.tsv.gz",
  "GRCh37_M11_O1.dedup_bcftools__CHR_nc_to_chr.nochr.aaf_scores.tsv.gz",
  "GRCh37_M12_C1.dedup_bcftools__CHR_nc_to_chr.nochr.aaf_scores.tsv.gz",
  "GRCh37_M12_O1.dedup_bcftools__CHR_nc_to_chr.nochr.aaf_scores.tsv.gz",
  "GRCh37_M13_C1.dedup_bcftools__CHR_nc_to_chr.nochr.aaf_scores.tsv.gz",
  "GRCh37_M13_O1.dedup_bcftools__CHR_nc_to_chr.nochr.aaf_scores.tsv.gz"
)

# Path to your phenocodes file (should contain 'filename' and 'description' columns)
phenocodes_path <- "phenocodes"

# Output composite plot filename
output_composite_plot_filename <- "PGS_Change_Dumbbell_Composite_PerPairTopPhenotypes_YBroken_XBroken_YReordered_CustomChunks.png"

# Trait selection parameters
PERCENTILE_THRESHOLD <- 99 # Filter by top N percentile of absolute differences (now per pair)
MAX_TRAITS_TO_PLOT <- 15    # Cap for readability (maximum phenotypes shown per donor)

# --- Y-Axis Breaking Configuration ---
# Default number of phenotypes to show per "chunk" if no custom setting for a donor
PHENOTYPES_PER_CHUNK <- 4

# CUSTOM CHUNK SIZES PER DONOR
# Specify the number of phenotypes for each chunk for specific donors.
# If a donor is not listed here, it will use the global PHENOTYPES_PER_CHUNK.
# The sum of chunk sizes for a donor should ideally not exceed MAX_TRAITS_TO_PLOT.
# Note: The script will adjust these if they try to plot more phenotypes than available
#        (up to MAX_TRAITS_TO_PLOT) for a given donor.
custom_chunk_sizes_per_donor <- list(
  "M8" = c(8, 4, 3),
  "M11" = c(2, 5, 7, 1),
  "M12" = c(1, 6, 6, 2),
  "M13" = c(1, 4, 4, 3, 2)
)

# --- CONFIGURABLE COLUMN NAMES IN YOUR INPUT TSV.GZ FILES ---
# !!! IMPORTANT: Adjust these to match the actual column names in your files !!!
PHENOTYPE_COL_NAME <- "match_file"
SCORE_COL_NAME <- "product_total"

# --- HELPER FUNCTIONS ---

#' Extracts the individual ID from a given file path.
#' Assumes ID is between the first underscore and the first dot.
extract_individual_id <- function(filepath) {
  filename <- basename(filepath)
  # 1. Get the part of the string after the first underscore
  after_first_underscore <- str_replace(filename, "^[^_]*_", "")
  # 2. From this new string, get the part before the first dot
  id <- str_replace(after_first_underscore, "\\..*$", "")
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
        # Construct the final Phenotype_Label using the cropped code and description
        Phenotype_Label = ifelse(
          !is.na(description),
          paste0(cropped_code_label, " \u2192 ", description), # \u2192 is the right arrow character
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

# --- Filter out specific phenotypes (Uncomment and customize if needed) ---
# message("\n--- Filtering out unwanted phenotypes ---")
# initial_pheno_count <- length(unique(all_raw_data_long$Phenotype_Label))
#
# # Define the pattern for terms to remove, separated by '|' for OR logic in regex
# terms_to_remove_pattern <- "Place of birth|television|Nap|Getting up|Z34|Alcohol|Time spent using computer|Hot drink|Lifetime number of sexual partners|fruit|Lamb|Worry|Cereal|PM: final answer|full time eduction|Cake intake|tanning|sunburn|physical activity|tiredness|shift work|Mother's age|body size at age 10|continuous-22501"
#
# all_raw_data_long <- all_raw_data_long %>%
#   filter(!str_detect(Phenotype_Label, terms_to_remove_pattern))
#
# filtered_pheno_count <- length(unique(all_raw_data_long$Phenotype_Label))
# message(paste0("  - Removed phenotypes containing: '", gsub("\\|", "', '", terms_to_remove_pattern), "'. Initial unique phenotypes: ", initial_pheno_count, ", After filter: ", filtered_pheno_count, "."))

message(paste0("  - Total raw data rows loaded: ", nrow(all_raw_data_long)))

# --- Generate Plots for Each Pair and Store in a List ---
message("\n--- Generating individual dumbbell plots for each donor pair ---")
plot_list <- list()

for (i in seq(1, length(input_files), by = 2)) {
  file_c_path <- input_files[i]
  file_o_path <- input_files[i + 1]
  
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
    plot_list[[paste0("plot_", donor_id)]] <- ggplot() +
      labs(title = paste0(donor_id, " (No Data)"),
           x = "PGS Score", y = "Phenotype") + # X and Y labels for placeholder
      theme_void() +
      theme(plot.title = element_text(hjust = 0.5, size = 12, face = "italic"),
            # Ensure titles are visible even on void theme
            axis.title.x = element_text(size = 10, face = "bold", color = "black"),
            axis.title.y = element_text(size = 10, face = "bold", color = "black"),
            plot.margin = unit(c(0.1, 0.25, 0.1, 0.25), "cm")) # Consistent margins
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
    message(paste0("    - No phenotypes met the ", PERCENTILE_THRESHOLD, "th percentile threshold for pair ", donor_id, ". Selecting top 1 trait as fallback."))
    selected_top_phenotypes_df <- head(current_pair_data %>% arrange(desc(absolute_difference)), 1)
  }
  
  message(paste0("    - Selected ", nrow(selected_top_phenotypes_df), " top phenotypes for pair ", donor_id, "."))
  
  # Order phenotypes by their minimum score for consistent Y-axis ordering
  phenotype_order_for_plot <- selected_top_phenotypes_df %>%
    mutate(min_score = pmin(Score_C, Score_O)) %>% # Calculate the minimum of the two scores
    arrange(min_score) %>% # Order by this minimum score
    pull(Phenotype_Label)
  
  selected_top_phenotypes_df$Phenotype_Label <- factor(selected_top_phenotypes_df$Phenotype_Label, levels = phenotype_order_for_plot)
  
  # --- Determine chunk sizes for the current donor ---
  current_donor_chunk_sizes <- NULL
  num_phenotypes_total <- nrow(selected_top_phenotypes_df) # Total phenotypes for this donor (up to MAX_TRAITS_TO_PLOT)
  
  if (donor_id %in% names(custom_chunk_sizes_per_donor)) {
    # If custom sizes are provided, use them but ensure they don't exceed available phenotypes
    custom_sizes_raw <- custom_chunk_sizes_per_donor[[donor_id]]
    
    adjusted_custom_sizes <- c()
    cumulative_phenos <- 0
    for (cs in custom_sizes_raw) {
      if (cs <= 0) {
        message(paste0("      Skipping non-positive custom chunk size (", cs, ") for donor ", donor_id, "."))
        next # Skip non-positive chunk sizes
      }
      
      phenos_to_take <- min(cs, num_phenotypes_total - cumulative_phenos)
      if (phenos_to_take > 0) {
        adjusted_custom_sizes <- c(adjusted_custom_sizes, phenos_to_take)
        cumulative_phenos <- cumulative_phenos + phenos_to_take
      }
      if (cumulative_phenos >= num_phenotypes_total) break # Stop if all phenotypes are assigned
    }
    current_donor_chunk_sizes <- adjusted_custom_sizes
    
    if (length(current_donor_chunk_sizes) == 0 && num_phenotypes_total > 0) {
      message(paste0("      Warning: Custom chunk sizes for donor ", donor_id, " are invalid or result in no chunks. Falling back to default."))
      # Fallback logic for default if custom failed
      num_full_chunks <- floor(num_phenotypes_total / PHENOTYPES_PER_CHUNK)
      remainder_chunk_size <- num_phenotypes_total %% PHENOTYPES_PER_CHUNK
      current_donor_chunk_sizes <- rep(PHENOTYPES_PER_CHUNK, num_full_chunks)
      if (remainder_chunk_size > 0) {
        current_donor_chunk_sizes <- c(current_donor_chunk_sizes, remainder_chunk_size)
      }
    } else {
      message(paste0("      Using adjusted custom chunk sizes for donor ", donor_id, ": ", paste(current_donor_chunk_sizes, collapse = ", ")))
    }
    
  } else {
    # Fallback to global PHENOTYPES_PER_CHUNK
    if (num_phenotypes_total > 0) {
      num_full_chunks <- floor(num_phenotypes_total / PHENOTYPES_PER_CHUNK)
      remainder_chunk_size <- num_phenotypes_total %% PHENOTYPES_PER_CHUNK
      
      current_donor_chunk_sizes <- rep(PHENOTYPES_PER_CHUNK, num_full_chunks)
      if (remainder_chunk_size > 0) {
        current_donor_chunk_sizes <- c(current_donor_chunk_sizes, remainder_chunk_size)
      }
    } else {
      current_donor_chunk_sizes <- c() # No phenotypes, no chunks
    }
    message(paste0("      Using default chunk size (", PHENOTYPES_PER_CHUNK, ") for donor ", donor_id, "."))
  }
  
  # --- Create sub-plots for Y-axis chunks using current_donor_chunk_sizes ---
  individual_donor_plots <- list()
  
  current_pheno_start_index <- 1
  chunk_counter <- 1
  
  # Loop through the determined chunk sizes
  for (chunk_size_for_this_iteration in current_donor_chunk_sizes) {
    if (current_pheno_start_index > num_phenotypes_total) {
      break # No more phenotypes to plot
    }
    
    # Calculate the end index for the current chunk, capped by total phenotypes
    chunk_end_index <- min(current_pheno_start_index + chunk_size_for_this_iteration - 1, num_phenotypes_total)
    
    # Select phenotypes for the current chunk
    chunk_phenotypes <- phenotype_order_for_plot[current_pheno_start_index:chunk_end_index]
    chunk_data <- selected_top_phenotypes_df %>%
      filter(Phenotype_Label %in% chunk_phenotypes)
    
    # Ensure chunk_data is not empty before proceeding with plotting for this chunk
    if (nrow(chunk_data) == 0) {
      message(paste0("        Warning: No data for chunk (", chunk_counter, ") of donor ", donor_id, ". Skipping chunk plot."))
      current_pheno_start_index <- chunk_end_index + 1 # Still advance index
      chunk_counter <- chunk_counter + 1
      next
    }
    
    # Prepare data for plotting points (long format) for this chunk
    chunk_plot_points_data <- chunk_data %>%
      pivot_longer(
        cols = c(Score_C, Score_O),
        names_to = "Score_Type",
        values_to = "Score_Value"
      ) %>%
      mutate(
        Dot_Color_Category = case_when(
          Score_Type == "Score_C" ~ "Centre",
          Score_Type == "Score_O" ~ "Outer"
        )
      )
    
    # Determine X-axis limits for this specific chunk plot based *only* on its data
    plot_min_score <- min(chunk_plot_points_data$Score_Value, na.rm = TRUE)
    plot_max_score <- max(chunk_plot_points_data$Score_Value, na.rm = TRUE)
    
    padding <- (plot_max_score - plot_min_score) * 0.1
    if (is.na(padding) || padding == 0) padding <- 1 # Ensure padding is not NA or zero, default to small value
    plot_min_score_padded <- plot_min_score - padding
    plot_max_score_padded <- plot_max_score + padding
    
    # Handle cases where min and max are the same (e.g., all scores are identical for this chunk)
    if (plot_min_score_padded == plot_max_score_padded) {
      plot_min_score_padded <- plot_min_score - 0.1
      plot_max_score_padded <- plot_max_score + 0.1
    }
    
    # DEBUG MESSAGE: Print the calculated X-axis range for this chunk
    message(paste0("        Chunk ", chunk_counter, " for donor ", donor_id, ": X-axis range [", round(plot_min_score_padded, 2), ", ", round(plot_max_score_padded, 2), "]"))
    
    p_chunk <- ggplot(chunk_data, aes(y = Phenotype_Label)) +
      geom_segment(aes(x = Score_C, xend = Score_O), linewidth = 0.8, color = "grey50") +
      geom_point(data = chunk_plot_points_data,
                 aes(x = Score_Value, color = Dot_Color_Category, shape = Dot_Color_Category),
                 size = 3, alpha=0.6) +
      labs(
        subtitle = if (chunk_counter == 1) paste0("Donor: ", donor_id) else NULL,
        x = "", # X-axis label on every chunk
        y = ""  # Y-axis label on every chunk
      ) +
      theme_minimal() +
      theme(
        plot.subtitle = element_text(hjust = 0.5, size = 12, face = "bold"),
        axis.title.x = element_text(size = 10, face = "bold", color = "black"), # Explicitly set X-axis title text properties
        axis.title.y = element_text(size = 10, face = "bold", color = "black"), # Explicitly set Y-axis title text properties
        axis.text.x = element_text(size = 8, color = "black", face = "bold"),
        axis.ticks.x = element_line(),
        axis.text.y = element_text(size = 6, color = "black", face = "bold"),
        plot.title = element_blank(), # Suppress individual chunk titles
        legend.position = "none", # Legends will be collected by patchwork globally
        panel.grid.major.y = element_line(color = "lightgray", linetype = "dotted", linewidth = 0.5),
        panel.grid.major.x = element_line(color = "lightgray", linetype = "dotted", linewidth = 0.5),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(0.1, 0.25, 0.1, 0.25), "cm") # Reduce margins between chunks
      ) +
      coord_cartesian(xlim = c(plot_min_score_padded, plot_max_score_padded)) +
      scale_color_manual(values = c("Centre" = "red", "Outer" = "blue")) +
      scale_shape_manual(values = c("Centre" = 16, "Outer" = 17))
    
    individual_donor_plots[[paste0("chunk_", chunk_counter)]] <- p_chunk
    
    # Move to the next chunk
    current_pheno_start_index <- chunk_end_index + 1
    chunk_counter <- chunk_counter + 1
  }
  
  if (length(individual_donor_plots) == 0) {
    message(paste0("    - No valid chunks created for donor pair ", donor_id, ". Adding empty placeholder."))
    plot_list[[paste0("plot_", donor_id)]] <- ggplot() +
      labs(title = paste0(donor_id, " (No Data)"), x = "", y = "") +
      theme_void() +
      theme(plot.title = element_text(hjust = 0.5, size = 12, face = "italic"),
            # Ensure titles are visible even on void theme for placeholder
            axis.title.x = element_text(size = 10, face = "bold", color = "black"),
            axis.title.y = element_text(size = 10, face = "bold", color = "black"),
            plot.margin = unit(c(0.1, 0.25, 0.1, 0.25), "cm")) # Consistent margins
  } else {
    # Combine chunks for the current donor vertically
    combined_donor_plot <- wrap_plots(individual_donor_plots, ncol = 1)
    
    # Add a common legend to the combined donor plot
    combined_donor_plot <- combined_donor_plot + plot_layout(guides = "collect") &
      scale_color_manual(name = "Sample Type",
                         values = c("Centre" = "red", "Outer" = "blue"),
                         labels = c(paste0("Centre (", sample_c_id, ")"), paste0("Outer (", sample_o_id, ")"))) &
      scale_shape_manual(name = "Sample Type",
                         values = c("Centre" = 16, "Outer" = 17),
                         labels = c(paste0("Centre (", sample_c_id, ")"), paste0("Outer (", sample_o_id, ")"))) &
      theme(legend.position = "bottom", legend.box = "horizontal", legend.direction = "horizontal")
    
    plot_list[[paste0("plot_", donor_id)]] <- combined_donor_plot
  }
}

# --- DEBUGGING: Check plot_list keys after loop ---
message("\n--- Debugging: Keys in plot_list after generation ---")
print(names(plot_list))
message(paste0("  - Total plots generated (by unique keys): ", length(plot_list)))
# --- END DEBUGGING ---

# --- Create Composite Plot ---
message("\n--- Creating composite plot ---")

n_plots <- length(plot_list)
n_cols_composite <- 2
n_rows_composite <- ceiling(n_plots / n_cols_composite)

composite_plot <- wrap_plots(plot_list, ncol = n_cols_composite) +
  plot_annotation(
    title = paste0("PGS Raw Score by Donor\nTop ", MAX_TRAITS_TO_PLOT, " Phenotypes, displaying ", PERCENTILE_THRESHOLD, "th percentile"),
    theme = theme(plot.title = element_text(hjust = 0.5, size = 11, face = "bold", margin = margin(b = 10)))
  ) & theme(legend.position = "bottom")

# Calculate dynamic height for the composite plot
total_height_per_plot_panel <- 0
for (plot_key in names(plot_list)) {
  donor_id <- gsub("plot_", "", plot_key)
  num_phenos_for_this_donor <- ifelse(donor_id %in% names(custom_chunk_sizes_per_donor),
                                      sum(custom_chunk_sizes_per_donor[[donor_id]]),
                                      MAX_TRAITS_TO_PLOT) # Approx if using default
  
  num_phenos_for_this_donor <- min(num_phenos_for_this_donor, MAX_TRAITS_TO_PLOT)
  
  # Rough estimate: 0.25 inch per phenotype, plus 2.0 inch for labels/margins/titles per donor plot
  total_height_per_plot_panel <- total_height_per_plot_panel + (num_phenos_for_this_donor * 0.25 + 2.0) # Increased base height slightly
}

composite_plot_height_in <- (total_height_per_plot_panel / n_cols_composite) + 3 # Add 3 for overall title/legend
composite_plot_width_in <- n_cols_composite * 8

# Cap max height/width to avoid excessively large plots
if (composite_plot_height_in > 60) composite_plot_height_in <- 60
if (composite_plot_width_in > 40) composite_plot_width_in <- 40

# --- SAVE THE COMPOSITE PLOT ---
ggsave(output_composite_plot_filename,
       plot = composite_plot,
       width = composite_plot_width_in,
       height = composite_plot_height_in,
       units = "in", dpi = 300, bg = "white")

message(paste0("\n--- Composite plot saved to '", output_composite_plot_filename, "' ---"))

message("--- Plotting complete ---")

composite_plot


ggsave(output_composite_plot_filename,
       plot = composite_plot,
       width = 18,  # <--- This is where the width is set
       height = 10, # <--- This is where the height is set
       units = "in", dpi = 300, bg = "white")
