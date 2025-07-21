Title: dumbbell_plot_composite_MC_%.R
Author: Dr. Alice M. Godden

## custom x axes for standardisation

# --- SCRIPT METADATA ---
# Project: PGS Score Change Visualization
# Description: Generates dumbbell plots to visualize PGS raw score changes for paired
#              samples (Control vs. Outcome). Features include dynamic Y-axis segmentation
#              based on percentage change range, X-axis adjusted per chunk, Y-axis reordering,
#              and explicit display of percentage change. Phenotype labels are now directly on the plot.

# Date: July 17, 2025 

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
setwd("~/PycharmProjects/polygenic/PGS_scores_test/SU/T0")


# --- CONFIGURATION PARAMETERS ---
# >>> IMPORTANT: Adjust these to match your actual file structure and desired behavior <<<

# Input files (raw score .tsv.gz files, assumed to be paired C, O, C, O...)
input_files <- c(
  #"GRCh37_M8_C1.dedup_bcftools__CHR_nc_to_chr.nochr.aaf_scores.tsv.gz",
  #"GRCh37_M8_O1.dedup_bcftools__CHR_nc_to_chr.nochr.aaf_scores.tsv.gz",
  "GRCh37_M11_C2._merged_dedup_bcftools__CHR_nc_to_chr.nochr.aaf_scores.tsv.gz",
  "GRCh37_M11_O1.dedup_bcftools__CHR_nc_to_chr.nochr.aaf_scores.tsv.gz",
  #"GRCh37_M12_C1.dedup_bcftools__CHR_nc_to_chr.nochr.aaf_scores.tsv.gz",
  #"GRCh37_M12_O1.dedup_bcftools__CHR_nc_to_chr.nochr.aaf_scores.tsv.gz"
  "GRCh37_M13_C1.dedup_bcftools__CHR_nc_to_chr.nochr.aaf_scores.tsv.gz",
  "GRCh37_M13_O1.dedup_bcftools__CHR_nc_to_chr.nochr.aaf_scores.tsv.gz"
)

# Path to your phenocodes file (should contain 'filename' and 'description' columns)
phenocodes_path <- "phenocodes"

# Output composite plot filename
output_composite_plot_filename <- "PGS_PercentageChange_Dumbbell_DynamicChunks_Composite.png"

# Trait selection parameters
PERCENTILE_THRESHOLD <- 99 # Filter by top N percentile of absolute differences (now per pair)
MAX_TRAITS_TO_PLOT <- 15    # Cap for readability (maximum phenotypes shown per donor)

# --- Dynamic Y-Axis Chunking Configuration (New) ---
# Maximum allowed range (difference between max and min percentage change) within a single plot chunk.
DYNAMIC_CHUNK_MAX_RANGE_PERCENT <- 75 # e.g., 75 percentage points. Chunks will split if the range exceeds this.

# Minimum and Maximum number of phenotypes per dynamic chunk.
# This helps prevent too many tiny chunks or single overly large chunks.
DYNAMIC_CHUNK_MIN_PHENO <- 1  # Aim for at least 1 phenotype per chunk
DYNAMIC_CHUNK_MAX_PHENO <- MAX_TRAITS_TO_PLOT # Max phenotypes per chunk set to MAX_TRAITS_TO_PLOT (15)

# --- X-Axis Standardization Parameters (Applied PER CHUNK) ---
# Define a minimum desired length for the X-axis of each chunk.
# If the actual data range within a chunk is smaller than this, the axis will extend to this minimum, centered.
MIN_X_AXIS_RANGE_PER_CHUNK <- 25 # e.g., 25% range minimum

# Define a fixed amount of additional breadth to add to the *effective* data range for each chunk,
# split equally on both sides.
ADDITIONAL_X_AXIS_BREADTH <- 20 # e.g., 20% fixed breadth

X_AXIS_PADDING_PERCENT <- 0.15 # Add 15% padding to the *effective* data range for each chunk

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
  donor_id <- str_extract(sample_c_id, "^M\\d+") # Extract donor ID (e.g., M8)
  
  message(paste0("  - Preparing plot for donor pair: ", donor_id, " (", sample_c_id, " vs ", sample_o_id, ")"))
  
  # Filter data for the current pair BEFORE pivoting
  current_pair_raw_filtered <- all_raw_data_long %>%
    filter(SampleID %in% c(sample_c_id, sample_o_id))
  
  # --- CRITICAL FIX: Add robustness check here for missing sample data ---
  if (!sample_c_id %in% unique(current_pair_raw_filtered$SampleID) ||
      !sample_o_id %in% unique(current_pair_raw_filtered$SampleID)) {
    missing_id <- if (!sample_c_id %in% unique(current_pair_raw_filtered$SampleID)) sample_c_id else sample_o_id
    message(paste0("    - Warning: Data for sample '", missing_id, "' is missing in the combined raw data. Skipping plot for donor ", donor_id, "."))
    plot_list[[paste0("plot_", donor_id)]] <- ggplot() +
      labs(title = paste0(donor_id, " (Data Missing for ", missing_id, ")"),
           x = "Percentage Change from Centre", y = "Phenotype") +
      theme_void() +
      theme(plot.title = element_text(hjust = 0.5, size = 12, face = "italic"),
            axis.title.x = element_text(size = 10, face = "bold", color = "black"),
            axis.title.y = element_text(size = 10, face = "bold", color = "black"),
            plot.margin = unit(c(0.1, 0.25, 0.1, 0.25), "cm"))
    next
  }
  # --- END CRITICAL FIX ---
  
  current_pair_data <- current_pair_raw_filtered %>%
    pivot_wider(names_from = SampleID, values_from = Score) %>%
    # Filter out rows where either C or O score is missing for the selected phenotypes
    # This filter will now correctly use the column names created by pivot_wider,
    # as the previous check ensures those columns exist (even if filled with NAs by pivot_wider)
    filter(!is.na(.data[[sample_c_id]]) & !is.na(.data[[sample_o_id]])) %>%
    mutate(
      Score_C = .data[[sample_c_id]],
      Score_O = .data[[sample_o_id]],
      absolute_difference = abs(Score_O - Score_C) # Keep for filtering
    )
  
  if (nrow(current_pair_data) == 0) {
    message(paste0("    - No complete data for donor pair ", donor_id, " after NA filtering. Skipping plot."))
    plot_list[[paste0("plot_", donor_id)]] <- ggplot() +
      labs(title = paste0(donor_id, " (No Data)"),
           x = "Percentage Change from Centre", y = "Phenotype") +
      theme_void() +
      theme(plot.title = element_text(hjust = 0.5, size = 12, face = "italic"),
            axis.title.x = element_text(size = 10, face = "bold", color = "black"),
            axis.title.y = element_text(size = 10, face = "bold", color = "black"),
            plot.margin = unit(c(0.1, 0.25, 0.1, 0.25), "cm"))
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
  
  # Order phenotypes by their percentage change for consistent Y-axis ordering within chunks
  # This is crucial for the dynamic chunking logic to work effectively.
  selected_top_phenotypes_df <- selected_top_phenotypes_df %>%
    arrange(Score_O_perc_change) %>% # Order by the percentage change of the outcome score
    mutate(Phenotype_Label = factor(Phenotype_Label, levels = unique(Phenotype_Label))) # Preserve order as factor levels
  
  # --- Dynamic Chunking Logic ---
  message(paste0("    - Applying dynamic chunking for donor ", donor_id, "..."))
  
  all_chunks_data <- list()
  current_chunk_phenos <- c()
  current_chunk_min_perc <- Inf
  current_chunk_max_perc <- -Inf
  
  if (nrow(selected_top_phenotypes_df) > 0) {
    for (row_idx in 1:nrow(selected_top_phenotypes_df)) {
      pheno <- selected_top_phenotypes_df[row_idx, ]
      pheno_label <- as.character(pheno$Phenotype_Label)
      perc_change <- pheno$Score_O_perc_change
      
      # Calculate potential new min/max if this phenotype is added to the current chunk
      potential_min_perc <- min(current_chunk_min_perc, perc_change)
      potential_max_perc <- max(current_chunk_max_perc, perc_change)
      
      start_new_chunk_flag <- FALSE
      
      if (length(current_chunk_phenos) == 0) { # Always start a new chunk for the first item
        start_new_chunk_flag <- TRUE
      } else if (length(current_chunk_phenos) >= DYNAMIC_CHUNK_MAX_PHENO) { # Current chunk is already "full" by count
        start_new_chunk_flag <- TRUE
      } else if ((potential_max_perc - potential_min_perc) > DYNAMIC_CHUNK_MAX_RANGE_PERCENT &&
                 length(current_chunk_phenos) >= DYNAMIC_CHUNK_MIN_PHENO) {
        # Adding this phenotype makes the range too large, and current chunk is already big enough
        start_new_chunk_flag <- TRUE
      }
      
      if (start_new_chunk_flag) {
        if (length(current_chunk_phenos) > 0) {
          # Save the previous chunk before starting a new one
          all_chunks_data[[length(all_chunks_data) + 1]] <- selected_top_phenotypes_df %>%
            filter(Phenotype_Label %in% current_chunk_phenos)
        }
        # Start new chunk with current phenotype
        current_chunk_phenos <- c(pheno_label)
        current_chunk_min_perc <- perc_change
        current_chunk_max_perc <- perc_change
      } else {
        # Add phenotype to current chunk
        current_chunk_phenos <- c(current_chunk_phenos, pheno_label)
        current_chunk_min_perc <- potential_min_perc
        current_chunk_max_perc <- potential_max_perc
      }
    }
    # Add the last chunk if not empty
    if (length(current_chunk_phenos) > 0) {
      all_chunks_data[[length(all_chunks_data) + 1]] <- selected_top_phenotypes_df %>%
        filter(Phenotype_Label %in% current_chunk_phenos)
    }
  }
  
  if (length(all_chunks_data) == 0) {
    message(paste0("    - No valid dynamic chunks created for donor pair ", donor_id, ". Adding empty placeholder."))
    plot_list[[paste0("plot_", donor_id)]] <- ggplot() +
      labs(title = paste0(donor_id, " (No Data)"), x = "Percentage Change from Centre (%)", y = "Phenotype") +
      theme_void() +
      theme(plot.title = element_text(hjust = 0.5, size = 12, face = "italic"),
            axis.title.x = element_text(size = 10, face = "bold", color = "black"),
            axis.title.y = element_text(size = 10, face = "bold", color = "black"),
            plot.margin = unit(c(0.1, 0.25, 0.1, 0.25), "cm"))
    next
  }
  
  # --- Create sub-plots for Y-axis chunks ---
  individual_donor_plots <- list()
  chunk_counter <- 1
  
  for (chunk_data in all_chunks_data) {
    # Re-order factor levels for the current chunk to match its data's sorted order
    chunk_data$Phenotype_Label <- factor(chunk_data$Phenotype_Label, levels = unique(chunk_data$Phenotype_Label))
    
    # Prepare data for plotting points (long format) for this chunk
    chunk_plot_points_data <- chunk_data %>%
      pivot_longer(
        cols = c(Score_C_perc_change, Score_O_perc_change), # Use percentage change columns
        names_to = "Score_Type",
        values_to = "Score_Value"
      ) %>%
      mutate(
        Dot_Color_Category = case_when(
          Score_Type == "Score_C_perc_change" ~ "Centre",
          Score_Type == "Score_O_perc_change" ~ "Outer"
        )
      )
    
    # Determine X-axis limits for this specific chunk plot based *only* on its data
    plot_min_perc <- min(chunk_plot_points_data$Score_Value, na.rm = TRUE)
    plot_max_perc <- max(chunk_plot_points_data$Score_Value, na.rm = TRUE)
    
    # Ensure 0 is included in the range for calculation
    plot_min_perc_inclusive_of_zero <- min(plot_min_perc, 0)
    plot_max_perc_inclusive_of_zero <- max(plot_max_perc, 0)
    
    actual_range <- plot_max_perc_inclusive_of_zero - plot_min_perc_inclusive_of_zero
    effective_range <- max(MIN_X_AXIS_RANGE_PER_CHUNK, actual_range) + ADDITIONAL_X_AXIS_BREADTH
    final_range <- effective_range * (1 + 2 * X_AXIS_PADDING_PERCENT)
    
    midpoint <- (plot_min_perc_inclusive_of_zero + plot_max_perc_inclusive_of_zero) / 2
    plot_min_perc_padded <- midpoint - (final_range / 2)
    plot_max_perc_padded <- midpoint + (final_range / 2)
    
    # Ensure padding is not zero if min/max are the same, or if data is very tight
    if (plot_min_perc_padded == plot_max_perc_padded) {
      plot_min_perc_padded <- plot_min_perc - 0.1
      plot_max_perc_padded <- plot_max_perc + 0.1
    }
    
    message(paste0("        Chunk ", chunk_counter, " for donor ", donor_id, ": X-axis PERCENTAGE range [", round(plot_min_perc_padded, 2), ", ", round(plot_max_perc_padded, 2), "] (Phenos: ", paste(chunk_data$Phenotype_Label, collapse = ", "), ")"))
    
    p_chunk <- ggplot(chunk_data, aes(y = Phenotype_Label)) +
      geom_segment(aes(x = Score_C_perc_change, xend = Score_O_perc_change), linewidth = 0.8, color = "grey50") +
      geom_point(data = chunk_plot_points_data,
                 aes(x = Score_Value, color = Dot_Color_Category, shape = Dot_Color_Category),
                 size = 3, alpha=0.9) +
      labs(
        subtitle = if (chunk_counter == 1) paste0("Donor: ", donor_id) else NULL,
        x = "Percentage Change from Centre (%)", # Consistent X-axis label
        y = "Phenotype" # Explicitly set Y-axis title here
      ) +
      theme_minimal() +
      theme(
        plot.subtitle = element_text(hjust = 0.5, size = 12, face = "bold"),
        axis.title.x = element_text(size = 10, face = "bold", color = "black"),
        axis.title.y = element_text(size = 10, face = "bold", color = "black"), # Ensure Y-axis title is visible
        axis.text.x = element_text(size = 8, color = "black", face = "bold"),
        axis.ticks.x = element_line(),
        axis.text.y = element_text(size = 10, color = "black", face = "bold"), # Y-axis label size to 6 for density
        plot.title = element_blank(), # Suppress individual chunk titles
        legend.position = "none", # Legends will be collected by patchwork globally
        legend.title = element_text(size = 9, face = "bold"),
        legend.text = element_text(size = 8, face="bold"),
        panel.grid.major.y = element_line(color = "lightgray", linetype = "dotted", linewidth = 0.5),
        panel.grid.major.x = element_line(color = "lightgray", linetype = "dotted", linewidth = 0.5),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(0.1, 0.25, 0.1, 0.25), "cm")
      ) +
      scale_color_manual(values = c("Centre" = "red", "Outer" = "blue")) +
      scale_shape_manual(values = c("Centre" = 16, "Outer" = 17))
    
    individual_donor_plots[[paste0("chunk_", chunk_counter)]] <- p_chunk
    chunk_counter <- chunk_counter + 1
  }
  
  if (length(individual_donor_plots) == 0) {
    message(paste0("    - No valid plots generated for donor pair ", donor_id, ". Adding empty placeholder."))
    plot_list[[paste0("plot_", donor_id)]] <- ggplot() +
      labs(title = paste0(donor_id, " (No Data)"), x = "Percentage Change from Centre (%)", y = "Phenotype") +
      theme_void() +
      theme(plot.title = element_text(hjust = 0.5, size = 12, face = "italic"),
            axis.title.x = element_text(size = 10, face = "bold", color = "black"),
            axis.title.y = element_text(size = 10, face = "bold", color = "black"),
            plot.margin = unit(c(0.1, 0.25, 0.1, 0.25), "cm"))
  } else {
    combined_donor_plot <- wrap_plots(individual_donor_plots, ncol = 1)
    
    combined_donor_plot <- combined_donor_plot + plot_layout(guides = "collect") &
      scale_color_manual(name = "Sample Type",
                         values = c("Centre" = "red", "Outer" = "blue"),
                         labels = c(paste0("Centre (", sample_c_id, ")"), paste0("Outer (", sample_o_id, ")"))) &
      scale_shape_manual(name = "Sample Type",
                         values = c("Centre" = 16, "Outer" = 17),
                         labels = c(paste0("Centre (", sample_c_id, ")"), paste0("Outer (", sample_o_id, ")"))) &
      theme(legend.position = "bottom", legend.box = "horizontal", legend.direction = "horizontal",
            legend.text = element_text(size = 8, face = "bold"))
    
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
    title = paste0("PGS Percentage Change by Donor (Relative to Centre)\nTop ", MAX_TRAITS_TO_PLOT, " Phenotypes (", PERCENTILE_THRESHOLD, "th percentile)"),
    theme = theme(plot.title = element_text(hjust = 0.5, size = 11, face = "bold", margin = margin(b = 10)))
  ) & theme(legend.position = "bottom")

# Calculate dynamic height for the composite plot
total_height_per_plot_panel <- 0
for (plot_key in names(plot_list)) {
  # Estimate chunk heights for the donor based on total phenotypes and dynamic chunking parameters
  donor_id <- gsub("plot_", "", plot_key)
  
  # Find the data for this donor to count actual phenotypes
  # This part of the height calculation still needs to assume MAX_TRAITS_TO_PLOT
  # if the actual data isn't easily accessible here from a higher scope.
  # For robustness, we will assume max possible if the specific data isn't trivially available.
  num_phenos_for_this_donor <- MAX_TRAITS_TO_PLOT # Assume max for height calculation to be safe
  
  # Calculate approximate number of chunks for this donor's plot based on DYNAMIC_CHUNK_MAX_PHENO
  num_chunks_for_this_donor <- ceiling(num_phenos_for_this_donor / DYNAMIC_CHUNK_MAX_PHENO)
  
  # Each phenotype row takes about 0.2 inches (estimated) + a bit for chunk spacing/titles
  # Plus a base height for each chunk
  estimated_height_per_chunk <- DYNAMIC_CHUNK_MAX_PHENO * 0.2 + 0.5 # 0.2in per phenotype, 0.5in for titles/padding
  total_height_per_plot_panel <- total_height_per_plot_panel + (estimated_height_per_chunk * num_chunks_for_this_donor)
}

# Add some buffer for overall title and legend
composite_plot_height_in <- (total_height_per_plot_panel / n_cols_composite) + 1.5 # 1.5 inches for overall title and legend
composite_plot_height_in <- max(composite_plot_height_in, 8) # Ensure a minimum height

message(paste0("  - Saving composite plot to '", output_composite_plot_filename, "'..."))
ggsave(output_composite_plot_filename,
       plot = composite_plot,
       width = n_cols_composite * 10, # Kept at 10 inches per column for ample width
       height = composite_plot_height_in,
       units = "in", dpi = 300, bg = "white")

message(paste0("\n--- Script finished. Composite plot saved to: ", output_composite_plot_filename, " ---"))
composite_plot
