#!/usr/bin/env Rscript

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)

args <- commandArgs(trailingOnly = TRUE)
host_status <- if (length(args) >= 2) args[2] else "not_used"

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

strip_extensions <- function(name) {
  repeat {
    new_name <- sub("\\.(fastq|fastp|fq|gz|meta|csv)$", "", name, perl = TRUE)
    if (new_name == name) break
    name <- new_name
  }
  return(name)
}

extract_sample_name <- function(filename) {
  # Remove all trailing extensions such as .fastq.gz or .fastp.fastq.gz
  name <- strip_extensions(filename)
  
  # Remove common processing suffixes
  name <- sub("_other$", "", name)
  name <- sub("_viral$", "", name)
  name <- sub("_classification_results$", "", name)
  name <- sub("_fixed$", "", name)
  name <- sub("_T1$", "", name)
  
  # Remove leftover '.fastp' strings in the middle of the name
  name <- sub("\\.fastp$", "", name)
  
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

human_depleted_reads_raw <- count_and_create_df("\\.fastq\\.gz$", dir = "nohuman")
human_depleted_present <- nrow(human_depleted_reads_raw) > 0
human_depleted_reads <- if (human_depleted_present) {
  human_depleted_reads_raw %>% rename(human_depleted = count)
} else {
  data.frame(sample = character(), human_depleted = numeric())
}

host_depleted_reads_raw <- if (dir.exists("nohost")) {
  count_and_create_df("\\.fastq\\.gz$", dir = "nohost")
} else {
  data.frame(sample = character(), count = numeric())
}
host_depleted_present <- nrow(host_depleted_reads_raw) > 0
host_depleted_reads <- if (host_depleted_present) {
  host_depleted_reads_raw %>% rename(host_depleted = count)
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

# Ensure expected columns exist even if their source directories were absent
expected_cols <- c("raw", "trimmed", "human_depleted", "host_depleted", "viral")
for (col in expected_cols) {
  if (!col %in% names(all_data)) {
    all_data[[col]] <- 0
  }
}

# Track which depletion stages actually ran
include_human <- host_status %in% c("human_only", "both") && human_depleted_present
include_other <- host_status %in% c("other_only", "both") && host_depleted_present

if (nrow(all_data) == 0) {
  cat("No data found. Check if files are present in the correct directories.\n")
} else {
  all_data <- all_data %>%
    mutate(across(everything(), ~replace_na(., 0)))
  
  if (!include_human) {
    all_data$human_depleted <- all_data$trimmed
  }
  if (!include_other) {
    all_data$host_depleted <- all_data$human_depleted
  }
  
  trimmed_reads_vals <- pmax(all_data$raw - all_data$trimmed, 0)
  human_reads_vals <- if (include_human) pmax(all_data$trimmed - all_data$human_depleted, 0) else rep(0, nrow(all_data))
  host_reads_vals  <- if (include_other) pmax(all_data$human_depleted - all_data$host_depleted, 0) else rep(0, nrow(all_data))
  non_viral_vals   <- if (include_other) {
    pmax(all_data$host_depleted - all_data$viral, 0)
  } else if (include_human) {
    pmax(all_data$human_depleted - all_data$viral, 0)
  } else {
    pmax(all_data$trimmed - all_data$viral, 0)
  }
  
  all_data <- all_data %>%
    mutate(
      trimmed_reads = trimmed_reads_vals,
      human_reads   = human_reads_vals,
      host_reads    = host_reads_vals,
      non_viral     = non_viral_vals
    )
  
  trimmed_reads_pct_vals <- round(ifelse(all_data$raw > 0, all_data$trimmed_reads / all_data$raw * 100, 0), 2)
  viral_pct_vals         <- round(ifelse(all_data$raw > 0, all_data$viral / all_data$raw * 100, 0), 2)
  non_viral_pct_vals     <- round(ifelse(all_data$raw > 0, all_data$non_viral / all_data$raw * 100, 0), 2)
  human_reads_pct_vals   <- if (include_human) round(ifelse(all_data$raw > 0, all_data$human_reads / all_data$raw * 100, 0), 2) else NULL
  host_reads_pct_vals    <- if (include_other) round(ifelse(all_data$raw > 0, all_data$host_reads / all_data$raw * 100, 0), 2) else NULL
  
  all_data <- all_data %>%
    mutate(
      trimmed_reads_pct = trimmed_reads_pct_vals,
      viral_pct         = viral_pct_vals,
      non_viral_pct     = non_viral_pct_vals
    )
  if (include_human) {
    all_data$human_reads_pct <- human_reads_pct_vals
  }
  if (include_other) {
    all_data$host_reads_pct <- host_reads_pct_vals
  }
  
  # Build ordered columns based on available host backgrounds
  column_order <- c("sample", "raw",
                    "trimmed_reads", "trimmed_reads_pct")
  if (include_human) {
    column_order <- c(column_order, "human_reads", "human_reads_pct")
  }
  if (include_other) {
    column_order <- c(column_order, "host_reads", "host_reads_pct")
  }
  column_order <- c(column_order, "viral", "viral_pct", "non_viral", "non_viral_pct")
  
  column_order <- column_order[column_order %in% names(all_data)]
  all_data <- all_data %>% select(all_of(column_order))
  
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
  
  available_categories <- unique(plot_data$category)
  category_order <- c("viral", "non_viral")
  if ("host_reads" %in% available_categories) {
    category_order <- c(category_order, "host_reads")
  }
  if ("human_reads" %in% available_categories) {
    category_order <- c(category_order, "human_reads")
  }
  category_order <- c(category_order, "trimmed_reads")
  category_order <- intersect(category_order, available_categories)
  plot_data$category <- factor(plot_data$category, levels = category_order)
  
  # Create color palette
  base_colors <- c("viral" = "#e78ac3",
                   "non_viral" = "#a6d854",
                   "host_reads" = "#8da0cb",
                   "human_reads" = "#fc8d62",
                   "trimmed_reads" = "#66c2a5")
  colors <- base_colors[names(base_colors) %in% category_order]
  
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
