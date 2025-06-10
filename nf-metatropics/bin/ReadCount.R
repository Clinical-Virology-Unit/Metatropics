#!/usr/bin/env Rscript

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)

# Set working directory to the 'read_count' folder
setwd("read_count")

# Print current working directory and list files
cat("Current working directory:", getwd(), "\n")
cat("Files in current directory:\n")
system("ls -R")

count_reads <- function(file_path) {
  tryCatch({
    con <- gzfile(file_path, "r")
    n <- 0
    while (length(readLines(con, n = 4)) > 0) {
      n <- n + 1
    }
    close(con)
    return(n)
  }, error = function(e) {
    warning(paste("Error reading file:", file_path, "-", e$message))
    return(0)
  })
}

extract_sample_name <- function(filename) {
  # Remove file extensions and common suffixes
  name <- sub("\\.[^.]+$", "", filename)  # Remove last extension
  name <- sub("_fixed$", "", name)        # Remove _fixed suffix
  name <- sub("_T1$", "", name)           # Remove _T1 suffix
  name <- sub("_other$", "", name)        # Remove _other suffix
  name <- sub("_viral$", "", name)        # Remove _viral suffix
  
  # Standardize sample names
  name <- gsub("H20_", "H2O_", name)      # Fix H20 to H2O
  name <- gsub("CRT57", "CRT-57", name)   # Fix CRT57 to CRT-57
  name <- gsub("UAC902", "UAC-902", name) # Fix UAC902 to UAC-902
  
  return(name)
}

count_and_create_df <- function(pattern, dir = ".") {
  cat("Searching for files matching pattern:", pattern, "in directory:", dir, "\n")
  
  # Make pattern more flexible for raw reads
  if (pattern == "^.*_T1\\.fastq\\.gz$") {
    pattern <- ".*_fixed\\.fastq\\.gz$"
  }
  
  files <- list.files(dir, pattern = pattern, full.names = TRUE)
  cat("Found", length(files), "files:\n")
  cat(paste(files, collapse = "\n"), "\n")
  
  if (length(files) == 0) {
    warning(paste("No files found matching pattern:", pattern, "in directory:", dir))
    return(data.frame(sample = character(), count = numeric()))
  }
  
  names <- basename(files)
  sample_names <- sapply(names, extract_sample_name)
  
  counts <- sapply(files, count_reads)
  df <- data.frame(sample = sample_names, count = counts)
  rownames(df) <- df$sample
  return(df)
}

raw_reads <- count_and_create_df("^.*_fixed\\.fastq\\.gz$")
if (nrow(raw_reads) > 0) raw_reads <- raw_reads %>% rename(raw = count)

trimmed_reads <- count_and_create_df("\\.fastp\\.fastq\\.gz$")
if (nrow(trimmed_reads) > 0) trimmed_reads <- trimmed_reads %>% rename(trimmed = count)

human_depleted_reads <- count_and_create_df("\\.fastq\\.gz$", dir = "nohuman")
if (nrow(human_depleted_reads) > 0) human_depleted_reads <- human_depleted_reads %>% rename(human_depleted = count)

host_depleted_reads <- if (dir.exists("nohost")) {
  df <- count_and_create_df("\\.fastq\\.gz$", dir = "nohost")
  if (nrow(df) > 0) df %>% rename(host_depleted = count) else df
} else {
  data.frame(sample = character(), host_depleted = numeric())
}

# Read viral reads from CSV file
viral_reads <- read.csv("viral_read_counts.csv", header = TRUE, stringsAsFactors = FALSE)
colnames(viral_reads) <- c("sample", "viral")
viral_reads$sample <- sapply(viral_reads$sample, extract_sample_name)
rownames(viral_reads) <- viral_reads$sample

all_data <- Reduce(function(x, y) full_join(x, y, by = "sample"),
                   list(raw_reads, trimmed_reads, human_depleted_reads, host_depleted_reads, viral_reads))

if (nrow(all_data) == 0) {
  cat("No data found. Check if files are present in the correct directories.\n")
} else {
  all_data <- all_data %>%
    mutate(across(everything(), ~replace_na(., 0))) %>%
    mutate(
      trimmed_reads = raw - trimmed,
      human_reads = trimmed - human_depleted
    )
  
  # Check if host depletion was performed
  host_depletion_performed <- dir.exists("nohost") && "host_depleted" %in% names(all_data) && sum(all_data$host_depleted) > 0
  
  if (host_depletion_performed) {
    all_data <- all_data %>%
      mutate(
        host_reads = human_depleted - host_depleted,
        non_viral = host_depleted - viral
      )
  } else {
    all_data <- all_data %>%
      mutate(
        non_viral = human_depleted - viral
      )
  }
  
  # Select columns based on whether host depletion was performed
  if (host_depletion_performed) {
    all_data <- all_data %>%
      select(sample, raw, trimmed_reads, human_reads, host_reads, viral, non_viral)
  } else {
    all_data <- all_data %>%
      select(sample, raw, trimmed_reads, human_reads, viral, non_viral)
  }
  
  # Calculate percentages
  all_data <- all_data %>%
    mutate(
      trimmed_reads_pct = round(trimmed_reads / raw * 100, 2),
      human_reads_pct = round(human_reads / raw * 100, 2),
      viral_pct = round(viral / raw * 100, 2),
      non_viral_pct = round(non_viral / raw * 100, 2)
    )
  
  # Add host_reads percentage only if host depletion was performed
  if (host_depletion_performed) {
    all_data <- all_data %>%
      mutate(host_reads_pct = round(host_reads / raw * 100, 2))
  }
  
  # Reorder columns
  if (host_depletion_performed) {
    all_data <- all_data %>%
      select(sample, raw, 
             trimmed_reads, trimmed_reads_pct, 
             human_reads, human_reads_pct, 
             host_reads, host_reads_pct, 
             viral, viral_pct, 
             non_viral, non_viral_pct)
  } else {
    all_data <- all_data %>%
      select(sample, raw, 
             trimmed_reads, trimmed_reads_pct, 
             human_reads, human_reads_pct, 
             viral, viral_pct, 
             non_viral, non_viral_pct)
  }
  
  rownames(all_data) <- all_data$sample
  all_data <- all_data %>% select(-sample)
  
  # Create a custom header with 'sample' as the first column name
  header <- c("sample", names(all_data))
  
  # Write the CSV file with the custom header
  write.table(rbind(header, cbind(rownames(all_data), all_data)), 
              file = "read_counts.csv", 
              sep = ",", 
              row.names = FALSE, 
              col.names = FALSE, 
              quote = FALSE)
  
  cat("Data processing completed. Results written to read_counts.csv\n")
  print(all_data)
  
  # Create stacked bar plot
  plot_data <- all_data %>%
    rownames_to_column("sample") %>%
    select(sample, ends_with("_pct")) %>%
    pivot_longer(cols = -sample, names_to = "category", values_to = "percentage") %>%
    mutate(category = sub("_pct$", "", category))
  
  # Define the new order of categories, accounting for optional host category
  if ("host_reads" %in% unique(plot_data$category)) {
    category_order <- c("viral", "non_viral", "host_reads", "human_reads", "trimmed_reads")
  } else {
    category_order <- c("viral", "non_viral", "human_reads", "trimmed_reads")
  }
  plot_data$category <- factor(plot_data$category, levels = category_order)
  
  # Create color palette
  colors <- c("viral" = "#e78ac3", "non_viral" = "#a6d854", 
              "host_reads" = "#8da0cb", "human_reads" = "#fc8d62", "trimmed_reads" = "#66c2a5")
  
  # Create the plot with improved aesthetics, borders around bars, and adjusted margins
  p <- ggplot(plot_data, aes(x = percentage, y = sample, fill = category)) +
    geom_bar(stat = "identity", color = "black", size = 0.25) +
    scale_fill_manual(values = colors, guide = guide_legend(reverse = TRUE)) +
    theme_bw() +
    theme(
      axis.text.y = element_text(angle = 0, hjust = 1),
      panel.border = element_rect(fill=NA, size=0.25),
      legend.position = "bottom",
      legend.title = element_blank(),
      plot.margin = margin(t = 10, r = 20, b = 10, l = 20, unit = "pt")
    ) +
    labs(title = "Read Distribution by Sample",
         x = "Percentage",
         y = NULL) +
    scale_x_continuous(labels = function(x) paste0(x, "%"), 
                       breaks = seq(0, 100, 25),
                       expand = c(0.01, 0)) +
    coord_cartesian(clip = "off")
  
  # Save the plot as PDF
  ggsave("read_distribution.pdf", plot = p, width = 10, height = 8)
  
  cat("Stacked bar plot saved as read_distribution.pdf\n")
}
