# polygenic risk score analysis for raw scores
setwd("~/PycharmProjects/polygenic/PGS_scores_test/SU/T0") # <<< ENSURE THIS PATH IS CORRECT

library(dplyr)      # For data manipulation (e.g., mutate, filter, group_by, summarise)
library(readr)      # For reading data (e.g., read_tsv)
library(ggplot2)    # For plotting
library(stringr)    # For string manipulation (e.g., str_replace, str_split)
library(tools)      # For file path manipulation (e.g., basename)
library(data.table) # For efficient data handling, especially if dealing with large files
library(patchwork)  # For combining multiple ggplot objects
library(scales)     # For log axis transformations (breaks_log, label_number)

# --- HELPER FUNCTION TO EXTRACT INDIVIDUAL ID FROM FILENAME ---
# MODIFIED: Extracts the segment between the first and second underscore.
extract_individual_id <- function(filepath) {
  filename <- basename(filepath)
  # Split the filename by underscore and take the second element.
  # Example: "GRCh37_D6T4a_merged..." -> Splits into ["GRCh37", "D6T4a", "merged...", ...]
  # We want the second element, which is "D6T4a".
  id <- str_split_i(filename, "_", 2) 
  return(id)
}

# --- MAIN DATA PROCESSING FUNCTION FOR AN INDIVIDUAL RAW SCORE FILE ---
# This function is adapted to read your _aaf_scores.tsv.gz files
process_raw_score_file <- function(filepath, verbose = TRUE) {
  # Reads an individual's raw score file.
  
  individual_id <- extract_individual_id(filepath)
  
  tryCatch({
    # read_tsv from readr package handles .gz files directly
    df_individual <- read_tsv(filepath, show_col_types = FALSE)
    
    if (nrow(df_individual) == 0) {
      if (verbose) message(paste("  - Warning: Raw score file", filepath, "is empty. Skipping."))
      return(NULL)
    }
    
    # *** CRITICAL LINE: Using "product_total" as identified from your file snippet ***
    score_col_name <- "product_total" 
    
    if (!(score_col_name %in% colnames(df_individual))) {
      if (verbose) {
        message(paste("  - Error: Missing expected score column ('", score_col_name, "') in raw score file", filepath, ". Skipping."))
        message(paste("    - Found columns:", paste(colnames(df_individual), collapse = ", ")))
      }
      return(NULL)
    }
    
    # Select only the score column and rename it to 'raw_score' for consistent internal processing
    df_individual_cleaned <- df_individual %>%
      select(raw_score = !!sym(score_col_name)) %>% # Use !!sym() to select by string name
      filter(!is.na(raw_score) & !is.infinite(raw_score) & !is.nan(raw_score)) # Remove non-finite values
    
    if (nrow(df_individual_cleaned) == 0) {
      if (verbose) message(paste("  - Warning: No valid raw score data in '", filepath, "' after cleaning. Skipping."))
      return(NULL)
    }
    
    df_individual_cleaned$individual_id <- individual_id
    # Return only the individual_id and the cleaned raw_score column
    return(df_individual_cleaned %>% select(individual_id, raw_score))
    
  }, error = function(e) {
    if (verbose) message(paste("  - An error occurred while reading raw score file", filepath, ":", e$message))
    return(NULL)
  })
}

# --- CONFIGURATION ---
# Specify all your raw score input files
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




# Output plot filename, updated to be more general for multiple plots
output_plot_filename <- "Raw_PGS_Score_Histograms_AllIndividuals_SU_RAW.png"

# --- PRE-ANALYSIS: Determine global x-axis range and global max frequency for common y-axis ---
message("\n--- Pre-analyzing files to determine global raw score range and global max frequency for common axes ---")
all_raw_scores <- numeric(0) # Will collect all raw scores from all files

if (length(input_files) == 0) {
  message("Error: No input files specified in the 'input_files' list. Exiting.")
  quit(save = "no", status = 1)
}

# Loop through files to gather all scores for global range determination
for (filepath in input_files) {
  processed_df <- process_raw_score_file(filepath, verbose = FALSE) # Read without verbose messages for pre-analysis
  if (!is.null(processed_df)) {
    all_raw_scores <- c(all_raw_scores, processed_df$raw_score)
  }
}

global_min_score <- min(all_raw_scores, na.rm = TRUE)
global_max_score <- max(all_raw_scores, na.rm = TRUE)

# Ensure a minimum range if all scores are identical or very close (avoids issues with binning)
if (global_max_score == global_min_score) {
  global_min_score <- global_min_score - 100
  global_max_score <- global_max_score + 100
  message("  - Adjusted global score range due to identical min/max values for plotting.")
}

message(paste0("  - Global Min Raw Score: ", round(global_min_score, 2)))
message(paste0("  - Global Max Raw Score: ", round(global_max_score, 2)))

# Define a common binwidth based on the global range. Aim for around 50 bins.
bin_width <- (global_max_score - global_min_score) / 50
if (bin_width == 0 || is.na(bin_width) || is.infinite(bin_width)) {
  bin_width <- 100 # Fallback for problematic ranges
  message(paste0("  - Warning: Calculated binwidth was problematic, falling back to: ", bin_width))
} else {
  if (bin_width < 1) bin_width <- 1 # Ensure bin_width isn't excessively small
}
message(paste0("  - Using histogram binwidth: ", round(bin_width, 2)))


# --- Calculate histogram counts for each individual with GLOBAL bin definitions ---
# This is crucial for a common Y-axis scale across plots
global_max_frequency <- 0 # To determine the upper limit for the Y-axis across all plots
all_binned_counts <- list() # To store binned counts for each individual file

for (filepath in input_files) {
  processed_df <- process_raw_score_file(filepath, verbose = FALSE)
  individual_id <- extract_individual_id(filepath)
  
  if (!is.null(processed_df) && nrow(processed_df) > 0) {
    # Calculate breaks based on global min/max score and bin_width
    breaks <- seq(floor(global_min_score / bin_width) * bin_width,
                  ceiling(global_max_score / bin_width) * bin_width,
                  by = bin_width)
    
    # Ensure breaks cover the full range of data, adding an extra bin if max is at the edge
    if (tail(breaks, 1) < global_max_score) {
      breaks <- c(breaks, breaks[length(breaks)] + bin_width)
    }
    
    # Determine bin lower boundaries
    bins_low <- breaks[-length(breaks)]
    
    # Initialize counts for all possible bins to zero
    counts <- rep(0, length(bins_low))
    
    # Efficiently assign scores to bins and count frequencies
    # findInterval returns the index of the interval a value falls into (1-based)
    bin_indices <- findInterval(processed_df$raw_score, bins_low, rightmost.closed = TRUE)
    
    # Aggregate counts for each bin index
    tabulated_counts <- table(bin_indices)
    
    # Update the 'counts' vector from the 'tabulated_counts', handling potential missing bins
    for (idx_str in names(tabulated_counts)) {
      idx <- as.numeric(idx_str)
      if (idx > 0 && idx <= length(counts)) { # Ensure index is valid
        counts[idx] <- counts[idx] + tabulated_counts[idx_str]
      }
    }
    
    # Create a dataframe of bin centers and their counts
    counts_df <- data.frame(
      bin_center = bins_low + bin_width / 2, # Midpoint of each bin
      count = counts
    )
    
    # Update global max frequency
    if (nrow(counts_df) > 0) {
      current_max_freq <- max(counts_df$count, na.rm = TRUE)
      if (current_max_freq > global_max_frequency) {
        global_max_frequency <- current_max_freq
      }
    }
    all_binned_counts[[individual_id]] <- counts_df # Store binned data for plotting
  } else {
    all_binned_counts[[individual_id]] <- data.frame(bin_center = numeric(0), count = numeric(0)) # Empty dataframe if no data
  }
}

# Determine the upper limit for the log scale y-axis based on global max frequency
if (global_max_frequency == 0) {
  y_axis_upper_limit_log <- 10 # Default if no data to avoid log(0) issues
  message("  - Warning: Global max frequency is 0. Setting log y-axis upper limit to 10.")
} else {
  max_freq_log_limit <- ceiling(log10(global_max_frequency + 1))
  y_axis_upper_limit_log <- 10^max_freq_log_limit
  if (y_axis_upper_limit_log <= 1) y_axis_upper_limit_log <- 10 # Ensure minimum range
}
y_axis_lower_limit_log <- 1 # Always start log Y at 1 for counts (10^0)

message(paste0("  - Global Max Frequency (for y-axis limit): ", global_max_frequency))
message(paste0("  - Log10 Y-axis Upper Limit: ", y_axis_upper_limit_log))

# --- PROCESS AND PLOT EACH INDIVIDUAL FILE ---
message("\n--- Processing and plotting individual files ---")

num_files <- length(input_files)

# Dynamic calculation of grid layout
ncols <- 2 # You can change this to 3 or more if you have many files and want a wider plot
nrows <- ceiling(num_files / ncols)

# Create a list to store ggplot objects
plot_list <- list()

for (i in 1:num_files) {
  filepath <- input_files[i]
  individual_id_for_plot <- extract_individual_id(filepath)
  
  # Retrieve pre-calculated binned counts for this file
  plot_data <- all_binned_counts[[individual_id_for_plot]]
  
  if (!is.null(plot_data) && nrow(plot_data) > 0) {
    # Replace 0 counts with NA so they are not drawn and don't cause issues with log scale
    plot_data_for_geom_bar <- plot_data %>% mutate(count_for_log = ifelse(count == 0, NA, count))
    
    p <- ggplot(plot_data_for_geom_bar, aes(x = bin_center, y = count_for_log)) +
      geom_col(width = bin_width, fill = "steelblue", color = "black") + # Histogram bars
      
      labs(title = paste0("", individual_id_for_plot),
           x = "Raw PGS Score",
           y = "Frequency (Log10 Scale)") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
            axis.title.x = element_text(size = 7, face = "bold", color = "black"),
            axis.title.y = element_text(size = 7, face = "bold", color = "black"),
            axis.text.x = element_text(size = 7, color = "black"),
            axis.text.y = element_text(size = 7, color = "black"),
            panel.grid.major.x = element_line(color = "lightgray", linetype = "dotted"),
            panel.grid.minor.x = element_blank(),
            panel.grid.major.y = element_line(color = "lightgray", linetype = "dotted"),
            panel.grid.minor.y = element_blank()) +
      
      # Set common x-axis limits for better comparison across plots
      coord_cartesian(xlim = c(global_min_score, global_max_score),
                      # Set common y-axis limits for log scale for consistency
                      ylim = c(y_axis_lower_limit_log, y_axis_upper_limit_log)) +
      scale_y_continuous(
        trans = 'log10', # Apply log10 transformation to the y-axis
        breaks = scales::breaks_log(n = 8), # Automatically determine nice log breaks
        labels = scales::label_number(accuracy = 1) # Format labels as numbers (1, 10, 100 etc.)
      )
    
    plot_list[[i]] <- p
  } else {
    # Placeholder plot if data processing failed for a file
    message(paste0("  - No valid raw scores or processed data found for '", filepath, "'. Skipping subplot or plotting empty placeholder."))
    plot_list[[i]] <- ggplot() +
      annotate("text", x=0.5, y=0.5, label="No data available or processing failed", size=5, color="red") +
      labs(title = paste0("Raw Score Histogram: ", individual_id_for_plot)) +
      theme_void()
  }
}

# Arrange plots into a grid
combined_plot <- wrap_plots(plot_list, ncol = ncols, nrow = nrows) +
  plot_annotation(
    title = paste0("Histograms of Raw PGS Scores"), # Main title for the combined plot
    theme = theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5, margin = margin(t = 5, b = 5)))
  ) &
  theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm")) # Adjust overall plot margins

# --- Save the composite plot ---
# Dimensions are now calculated dynamically based on the number of rows and columns
subplot_width_in <- 7        # Width of a single subplot in inches (adjust as needed)
subplot_height_in <- 4.5     # Height of a single subplot in inches (adjust as needed)

final_plot_width <- ncols * subplot_width_in   # Total width of the combined plot
final_plot_height <- nrows * subplot_height_in # Total height of the combined plot

ggsave(output_plot_filename,
       plot = combined_plot,
       width = final_plot_width,
       height = final_plot_height,
       units = "in", dpi = 600, bg = "white") # High resolution for quality

message(paste0("\n--- Composite raw score histogram plot saved to '", output_plot_filename, "' ---"))

# To display the plot (this will usually happen automatically in RStudio's Plots pane)
print(combined_plot)

message("--- Plotting complete ---")