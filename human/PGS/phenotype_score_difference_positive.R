# polygenic risk score analysis
setwd("~/PycharmProjects/polygenic/PGS_scores_test/SU/T0")



library(dplyr)      # For data manipulation (e.g., mutate, filter, group_by, summarise)
library(readr)      # For reading data (e.g., read_tsv)
library(ggplot2)    # For plotting
library(stringr)    # For string manipulation (e.g., str_replace, str_split)
library(tools)      # For file path manipulation (e.g., basename)
library(data.table) # For efficient data handling, especially if dealing with large files
library(patchwork)  # For combining multiple ggplot objects

# --- HELPER FUNCTION TO CLEAN TRAIT LABELS ---
clean_label <- function(label) {
  # Cleans the filename to extract a more readable trait label.
  label_str <- as.character(label)
  cleaned_label <- str_replace_all(label_str, "_filtered.tsv.bgz", "")
  cleaned_label <- str_split(cleaned_label, "-both_sexes", simplify = TRUE)[, 1]
  
  # MODIFICATION: Refined logic to restrict label to maximum of second underscore
  # This tries to get the first two parts if available, otherwise just keeps what's there
  temp_parts <- str_split(cleaned_label, "_", simplify = TRUE)
  if (ncol(temp_parts) >= 2) { # If there's at least one underscore (so at least 2 parts)
    cleaned_label <- paste(temp_parts[,1], temp_parts[,2], sep = "_")
  } else { # If there's 0 or 1 underscore, keep the original cleaned_label (first part)
    cleaned_label <- temp_parts[,1]
  }
  
  return(cleaned_label)
}


# --- MAIN DATA PROCESSING FUNCTION FOR AN INDIVIDUAL FILE ---
process_individual_file <- function(filepath, phenocodes_df, verbose = TRUE) {
  # Reads an individual's score difference file, cleans labels, and merges phenocodes.
  # ASSUMPTION: File contains match_file and score_difference for one individual.
  
  individual_id <- tools::file_path_sans_ext(basename(filepath)) # Get base name and remove extension
  individual_id <- str_replace(individual_id, "PGSscores_difference_", "") # Remove prefix
  
  tryCatch({
    df_individual <- read_tsv(filepath, show_col_types = FALSE) # read_tsv from readr package
    
    if (nrow(df_individual) == 0) {
      if (verbose) message(paste("  - Warning: Individual file", filepath, "is empty. Skipping."))
      return(NULL)
    }
    
    required_cols <- c('match_file', 'score_difference')
    if (!all(required_cols %in% colnames(df_individual))) {
      if (verbose) {
        message(paste("  - Error: Missing required columns (", paste(required_cols, collapse = ", "), ") in individual file", filepath, ". Skipping."))
        message(paste("    - Found columns:", paste(colnames(df_individual), collapse = ", ")))
      }
      return(NULL)
    }
    
    df_individual <- df_individual %>%
      mutate(label = sapply(match_file, clean_label), # Apply clean_label function
             filename_clean = str_replace_all(match_file, "_filtered.tsv.bgz", ".tsv.bgz"))
    
    # Merge phenocodes
    merged_df <- left_join(df_individual, phenocodes_df, by = c("filename_clean" = "filename"))
    
    # Only take the first part of the phenocode before the arrow, and then apply the short label logic
    merged_df <- merged_df %>%
      mutate(full_label_with_desc = ifelse(!is.na(description), paste0(label, " \u2192 ", description), label),
             label_short = str_split_i(full_label_with_desc, " \u2192 ", 1)) # Extract part before arrow
    
    # --- DEBUGGING FOR MC_M8_C-O SPECIFICALLY (removed some verbose for clarity here, but keep in your code) ---
    if (individual_id == "MC_M8_C-O" && verbose) {
      message("--- Debugging MC_M8_C-O label processing (before cleaning non-finite values) ---")
      message("  - Before score_difference conversion and cleaning, head(merged_df):")
      print(head(merged_df %>% select(match_file, label, full_label_with_desc, label_short, score_difference)))
      message("  - Check original score_difference values for MC_M8_C-O (summary):")
      print(summary(merged_df$score_difference))
      message("  - Check original score_difference values for MC_M8_C-O (top 5 values):")
      print(head(merged_df %>% arrange(desc(score_difference)) %>% select(label_short, score_difference)))
    }
    # --- END DEBUGGING ---
    
    # Convert score_difference to numeric, coercing errors to NA
    merged_df$score_difference <- as.numeric(merged_df$score_difference)
    
    # NEW: Drop rows where score_difference is Inf, -Inf, or NaN, in addition to NA
    # This specifically addresses the 'Inf' issue.
    merged_df_cleaned <- merged_df %>%
      filter(!is.na(score_difference) & !is.infinite(score_difference) & !is.nan(score_difference) &
               !is.na(label_short) & label_short != "")
    
    # --- DEBUGGING FOR MC_M8_C-O SPECIFICALLY AFTER FILTERING ---
    if (individual_id == "MC_M8_C-O" && verbose) {
      message("  - After score_difference conversion and cleaning (and Inf/NaN removal), head(merged_df_cleaned):")
      print(head(merged_df_cleaned %>% select(match_file, label, full_label_with_desc, label_short, score_difference)))
      message(paste("  - Rows in merged_df_cleaned:", nrow(merged_df_cleaned)))
      message("  - Summary of score_difference in merged_df_cleaned (after Inf/NaN removal):")
      print(summary(merged_df_cleaned$score_difference))
      message("  - Top 5 positive score_difference values after cleaning and Inf/NaN removal:")
      print(head(merged_df_cleaned %>% filter(score_difference > 0) %>% arrange(desc(score_difference)) %>% select(label_short, score_difference)))
    }
    # --- END DEBUGGING ---
    
    if (nrow(merged_df_cleaned) == 0) {
      if (verbose) message(paste("  - Warning: No valid score_difference data in '", filepath, "' after cleaning. Skipping."))
      return(NULL)
    }
    
    # Check for duplicate labels and take mean (now based on label_short)
    if (any(duplicated(merged_df_cleaned$label_short))) {
      if (verbose) message(paste("    - Warning: Duplicate trait labels found in", filepath, ". Taking mean score_difference for duplicates."))
      merged_df_cleaned <- merged_df_cleaned %>%
        group_by(label_short) %>%
        summarise(score_difference = mean(score_difference, na.rm = TRUE)) %>%
        ungroup() # Ungroup after summarising
    }
    
    # --- DEBUGGING FOR MC_M8_C-O SPECIFICALLY AFTER DUPLICATE HANDLING ---
    if (individual_id == "MC_M8_C-O" && verbose) {
      message("  - After duplicate handling, head(merged_df_cleaned):")
      print(head(merged_df_cleaned %>% select(label_short, score_difference)))
      message(paste("  - Rows in merged_df_cleaned after duplicate handling:", nrow(merged_df_cleaned)))
      message("--- End Debugging MC_M8_C-O label processing ---")
    }
    # --- END DEBUGGING ---
    
    if (any(duplicated(merged_df_cleaned$label_short))) { # Re-check for uniqueness after aggregation
      if (verbose) message(paste("    - Critical Error: Trait labels are still not unique in", filepath, "after mean aggregation. This should not happen. Skipping."))
      return(NULL)
    }
    
    merged_df_cleaned$individual_id <- individual_id
    # Return label_short for plotting
    return(merged_df_cleaned %>% select(individual_id, label_short, score_difference))
    
  }, error = function(e) {
    if (verbose) message(paste("  - An error occurred while reading individual file", filepath, ":", e$message))
    return(NULL)
  })
}

# ... (rest of the script remains the same from CONFIGURATION onwards) ...
    required_cols <- c('match_file', 'score_difference')
    if (!all(required_cols %in% colnames(df_individual))) {
      if (verbose) {
        message(paste("  - Error: Missing required columns (", paste(required_cols, collapse = ", "), ") in individual file", filepath, ". Skipping."))
        message(paste("    - Found columns:", paste(colnames(df_individual), collapse = ", ")))
      }
      return(NULL)
    }
    
    df_individual <- df_individual %>%
      mutate(label = sapply(match_file, clean_label), # Apply clean_label function
             filename_clean = str_replace_all(match_file, "_filtered.tsv.bgz", ".tsv.bgz"))
    
    # Merge phenocodes
    merged_df <- left_join(df_individual, phenocodes_df, by = c("filename_clean" = "filename"))
    
    # Only take the first part of the phenocode before the arrow, and then apply the short label logic
    merged_df <- merged_df %>%
      mutate(full_label_with_desc = ifelse(!is.na(description), paste0(label, " \u2192 ", description), label),
             label_short = str_split_i(full_label_with_desc, " \u2192 ", 1)) # Extract part before arrow
    
    # --- DEBUGGING FOR MC_M8_C-O SPECIFICALLY ---
    if (individual_id == "MC_M8_C-O" && verbose) {
      message("--- Debugging MC_M8_C-O label processing ---")
      message("  - Before score_difference conversion and cleaning, head(merged_df):")
      print(head(merged_df %>% select(match_file, label, full_label_with_desc, label_short, score_difference)))
      message("  - Check for NAs/empty strings in label_short:")
      message(paste("    - NAs in label_short:", sum(is.na(merged_df$label_short))))
      message(paste("    - Empty strings in label_short:", sum(merged_df$label_short == "", na.rm = TRUE)))
      message("  - Check original score_difference values for MC_M8_C-O:")
      print(summary(merged_df$score_difference))
    }
    # --- END DEBUGGING ---
    
    # Convert score_difference to numeric, coercing errors to NA
    merged_df$score_difference <- as.numeric(merged_df$score_difference)
    
    # Drop rows where score_difference or label is NA, or label_short is empty
    merged_df_cleaned <- merged_df %>%
      filter(!is.na(score_difference) & !is.na(label_short) & label_short != "")
    
    # --- DEBUGGING FOR MC_M8_C-O SPECIFICALLY AFTER FILTERING ---
    if (individual_id == "MC_M8_C-O" && verbose) {
      message("  - After score_difference conversion and cleaning, head(merged_df_cleaned):")
      print(head(merged_df_cleaned %>% select(match_file, label, full_label_with_desc, label_short, score_difference)))
      message(paste("  - Rows in merged_df_cleaned:", nrow(merged_df_cleaned)))
      message("  - Summary of score_difference in merged_df_cleaned:")
      print(summary(merged_df_cleaned$score_difference))
    }
    # --- END DEBUGGING ---
    
    
    if (nrow(merged_df_cleaned) == 0) {
      if (verbose) message(paste("  - Warning: No valid score_difference data in '", filepath, "' after cleaning. Skipping."))
      return(NULL)
    }
    
    # Check for duplicate labels and take mean (now based on label_short)
    if (any(duplicated(merged_df_cleaned$label_short))) {
      if (verbose) message(paste("    - Warning: Duplicate trait labels found in", filepath, ". Taking mean score_difference for duplicates."))
      merged_df_cleaned <- merged_df_cleaned %>%
        group_by(label_short) %>%
        summarise(score_difference = mean(score_difference, na.rm = TRUE)) %>%
        ungroup() # Ungroup after summarising
    }
    
    # --- DEBUGGING FOR MC_M8_C-O SPECIFICALLY AFTER DUPLICATE HANDLING ---
    if (individual_id == "MC_M8_C-O" && verbose) {
      message("  - After duplicate handling, head(merged_df_cleaned):")
      print(head(merged_df_cleaned %>% select(label_short, score_difference)))
      message(paste("  - Rows in merged_df_cleaned after duplicate handling:", nrow(merged_df_cleaned)))
      message("--- End Debugging MC_M8_C-O label processing ---")
    }
    # --- END DEBUGGING ---
    
    if (any(duplicated(merged_df_cleaned$label_short))) { # Re-check for uniqueness after aggregation
      if (verbose) message(paste("    - Critical Error: Trait labels are still not unique in", filepath, "after mean aggregation. This should not happen. Skipping."))
      return(NULL)
    }
    
    merged_df_cleaned$individual_id <- individual_id
    # Return label_short for plotting
    return(merged_df_cleaned %>% select(individual_id, label_short, score_difference))
    
  }, error = function(e) {
    if (verbose) message(paste("  - An error occurred while reading individual file", filepath, ":", e$message))
    return(NULL)
  })
}

# --- CONFIGURATION ---
#input_files <- c( # MC files
#  "PGSscores_difference_MC_M8_C-O.txt",
#  "PGSscores_difference_MC_M11_C-O.txt",
#  "PGSscores_difference_MC_M12_C-O.txt",
#  "PGSscores_difference_MC_M13_C-O.txt"
#)

# SwimUp files
input_files <- c(
"PGSscores_difference_D1T4-D1T0.txt",
"PGSscores_difference_D2T4-D2T0.txt",
"PGSscores_difference_D4T4-D4T0.txt",
"PGSscores_difference_D6aT4-D6T0a.txt",
"PGSscores_difference_D6bT4-D6bT0.txt"
)
phenocodes_path <- "phenocodes"
output_plot_filename <- "PGSscore_positive_differences__annotated_top40_SU.png" # Updated filename
MAX_TRAITS_TO_DISPLAY <- 40 # Changed back to 40

# --- LOAD PHENOCODES ---
message(paste0("\n--- Loading phenocodes from: ", phenocodes_path, " ---"))
tryCatch({
  phenos <- read_tsv(phenocodes_path, show_col_types = FALSE)
  message(paste0("  - Successfully loaded phenocodes. Dimensions: ", nrow(phenos), " rows, ", ncol(phenos), " columns"))
  if (!("filename" %in% colnames(phenos)) || !("description" %in% colnames(phenos))) {
    message(paste0("  - Error: 'filename' or 'description' column missing in '", phenocodes_path, "'. Please check phenocodes file."))
    quit(save = "no", status = 1) # Exit R script
  }
}, error = function(e) {
  message(paste0("Error loading phenocodes: ", e$message))
  quit(save = "no", status = 1) # Exit R script
})

# --- PRE-ANALYSIS: DETERMINE MAX NUMBER OF *DISPLAYED* PHENOTYPES FOR X-AXIS LIMIT AND COUNT TRAITS ---
message("\n--- Pre-analyzing files to count traits above 0 and determine max X-axis size (limited by MAX_TRAITS_TO_DISPLAY) ---")
actual_max_phenotypes_after_limit <- 0
trait_counts_per_donor <- list()

if (length(input_files) == 0) {
  message("Error: No input files specified in the 'input_files' list. Exiting.")
  quit(save = "no", status = 1)
}

for (filepath in input_files) {
  processed_df <- process_individual_file(filepath, phenos, verbose = FALSE) # Set verbose to FALSE for pre-analysis
  individual_id <- tools::file_path_sans_ext(basename(filepath))
  individual_id <- str_replace(individual_id, "PGSscores_difference_", "")
  
  if (!is.null(processed_df)) {
    positive_diff_df <- processed_df %>% filter(score_difference > 0)
    full_count_positive <- nrow(positive_diff_df)
    trait_counts_per_donor[[individual_id]] <- full_count_positive
    
    num_to_display <- min(full_count_positive, MAX_TRAITS_TO_DISPLAY)
    if (num_to_display > actual_max_phenotypes_after_limit) {
      actual_max_phenotypes_after_limit <- num_to_display
    }
  } else {
    trait_counts_per_donor[[individual_id]] <- 0
    message(paste0("  - Warning: Could not process '", filepath, "' during pre-analysis. Assuming 0 positive traits."))
  }
}

message(paste0("  - Max actual phenotypes to display on any single graph (after limit of ", MAX_TRAITS_TO_DISPLAY, "): ", actual_max_phenotypes_after_limit))

message("\n--- Summary of ALL traits above 0 per donor (before display limit) ---")
for (donor in names(trait_counts_per_donor)) {
  message(paste0("  ", donor, ": ", trait_counts_per_donor[[donor]]))
}

# Adjust figure height based on the maximum number of traits to display
base_fig_height <- 2
height_per_trait <- 0.35
figure_height <- max(base_fig_height, actual_max_phenotypes_after_limit * height_per_trait)
if (figure_height < 4) figure_height <- 4 # Ensure minimum height

# --- PROCESS AND PLOT EACH INDIVIDUAL FILE ---
message("\n--- Processing and plotting individual files ---")

num_files <- length(input_files)
ncols <- 2
nrows <- ceiling(num_files / ncols)

# Create a list to store plots
plot_list <- list()

for (i in 1:num_files) {
  filepath <- input_files[i]
  if (file.exists(filepath)) {
    # Set verbose to TRUE only for the actual plotting loop for specific debugging messages
    processed_df <- process_individual_file(filepath, phenos, verbose = TRUE)
    
    if (!is.null(processed_df)) {
      positive_diff_df <- processed_df %>% filter(score_difference > 0)
      individual_id_for_plot <- if (nrow(positive_diff_df) > 0) positive_diff_df$individual_id[1] else str_replace(tools::file_path_sans_ext(basename(filepath)), "PGSscores_difference_", "")
      num_positive_phenotypes_for_display <- nrow(positive_diff_df) # Total positive count
      
      if (nrow(positive_diff_df) > 0) {
        positive_diff_df_to_plot <- positive_diff_df %>%
          arrange(desc(score_difference)) %>% # Still sort by score_difference
          head(MAX_TRAITS_TO_DISPLAY) # Limit for plotting
        
        # Ensure labels are factors in the desired order
        positive_diff_df_to_plot$label_short <- factor(positive_diff_df_to_plot$label_short, levels = positive_diff_df_to_plot$label_short)
        
        # Dynamic Y-AXIS: Determine max score difference for this specific plot
        max_score_diff_for_plot <- max(positive_diff_df_to_plot$score_difference, na.rm = TRUE)
        # Add a small buffer to the max for visual comfort
        y_axis_upper_limit <- max_score_diff_for_plot * 1.1 # Increased buffer slightly for annotations
        
        p <- ggplot(positive_diff_df_to_plot, aes(x = label_short, y = score_difference)) + # Use label_short for x-axis
          geom_bar(stat = "identity", fill = "steelblue") +
          coord_cartesian(ylim = c(0, y_axis_upper_limit)) + # Dynamic Y-axis limits
          labs(title = paste0("Positive PGS Differences: ", individual_id_for_plot),
               x = "", # No label for x-axis as it's categorical traits
               y = "PGS Score Difference") + # Y-axis label for score difference
          theme_minimal() +
          theme(plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
                axis.title.y = element_text(size = 10, face = "bold", color="black"), # Y-axis title
                axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 6, color="black", face="bold"), # Smaller font size (e.g., 6)
                axis.ticks.x = element_line(), # Keep x-axis ticks visible
                panel.grid.major.x = element_blank(), # Remove vertical grid lines
                panel.grid.minor.x = element_blank(),
                axis.text.y = element_text(size = 8, color="black", face="bold")) +
          # MODIFICATION: Add count annotations back
          annotate("text",
                   x = Inf, y = Inf, # Position at top-right corner of plot
                   label = paste0("Phenotypes with PGS >0: ", num_positive_phenotypes_for_display),
                   hjust = 1.05, vjust = 1.2, # Adjust justification to move slightly inside
                   size = 3.5, color = "dimgray",
                   fontface = "plain")
        
        #if (num_positive_phenotypes_for_display > MAX_TRAITS_TO_DISPLAY) {
         # p <- p + annotate("text",
          #                  x = Inf, y = Inf, # Position at top-right
           #                 label = paste0("Showing top ", nrow(positive_diff_df_to_plot)),
            #                hjust = 1.05, vjust = 2.4, # Adjust justification (lower than first annotation)
             #               size = 2.5, color = "red",
              #              fontface = "plain")
        #}
        
        plot_list[[i]] <- p
      } else {
        message(paste0("  - No positive score differences found for '", filepath, "'. Skipping subplot."))
        plot_list[[i]] <- ggplot() + theme_void() # Create an empty plot to occupy the space
      }
    } else {
      message(paste0("  - Could not process '", filepath, "'. Skipping subplot."))
      plot_list[[i]] <- ggplot() + theme_void() # Create an empty plot to occupy the space
    }
  } else {
    message(paste0("Warning: Input file not found: '", filepath, "'. Skipping subplot."))
    plot_list[[i]] <- ggplot() + theme_void() # Create an empty plot to occupy the space
  }
}

# Arrange plots into a grid
combined_plot <- wrap_plots(plot_list, ncol = ncols, nrow = nrows) +
  plot_annotation(
    title = paste0("Positive PGS Score Differences per Individual (Showing Top ", MAX_TRAITS_TO_DISPLAY, " Traits)"),
    theme = theme(plot.title = element_text(size = 11, face = "bold", hjust = 0.5, margin = margin(t = 5, b = 5)))
  ) &
  theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm")) # Adjust overall plot margins

# Save the composite plot
ggsave(output_plot_filename, plot = combined_plot, width = ncols * 7, height = nrows * figure_height, units = "in", dpi = 300, bg = "white")
message(paste0("\n--- Composite plot saved to '", output_plot_filename, "' ---"))

# To display the plot (this will usually happen automatically in RStudio's Plots pane)
print(combined_plot)

message("--- Plotting complete ---")
# Save the composite plot
ggsave(output_plot_filename, plot = combined_plot, width = ncols * 7, height = nrows * figure_height, units = "in", dpi = 600, bg = "white")
message(paste0("\n--- Composite plot saved to '", output_plot_filename, "' ---"))

# To display the plot (this will usually happen automatically in RStudio's Plots pane)
print(combined_plot)
