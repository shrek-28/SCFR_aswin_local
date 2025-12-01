#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(ggplot2)
library(readr)
library(dplyr)

# Function to convert bp to bp/kb/Mb with 2 decimals
fmt_bp <- function(x) {
  if (x >= 1e6) {
    paste0(format(round(x / 1e6, 2), nsmall = 2), " Mb")
  } else if (x >= 1e3) {
    paste0(format(round(x / 1e3, 2), nsmall = 2), " Kb")
  } else {
    paste0(format(round(x, 2), nsmall = 2), " bp")
  }
}

# Load data
df <- read_tsv(args[1])

# Compute stats
stats_df <- df %>%
  group_by(species) %>%
  summarise(
    mean = mean(length),
    median = median(length),
    mode = as.numeric(names(sort(table(length), decreasing = TRUE))[1]),
    q1 = quantile(length, 0.25),
    q3 = quantile(length, 0.75),
    sd = sd(length)
  )

# Apply human-readable formatting
stats_df <- stats_df %>%
  mutate(
    label = paste0(
      "Mean: ", fmt_bp(mean), "\n",
      "Median: ", fmt_bp(median), "\n",
      "Mode: ", fmt_bp(mode), "\n",
      "Q1: ", fmt_bp(q1), "\n",
      "Q3: ", fmt_bp(q3), "\n",
      "SD: ", fmt_bp(sd)
    )
  )

# Boxplot with log10 y-scale
p <- ggplot(df, aes(x = species, y = length, fill = species)) +
  geom_boxplot(outlier.shape = 21, alpha = 0.85) +
  scale_y_log10() +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  labs(
    title = "Gene Desert Length Distribution Across Species",
    x = "Species",
    y = "Desert Length (log10 scale)"
  ) +
  geom_text(
    data = stats_df,
    aes(
      x = species,
      y = max(df$length) * 1.5,
      label = label
    ),
    size = 3.5,
    hjust = 0.5
  )

# Save as high-resolution PDF
ggsave("gene_desert_boxplot.pdf", p, width = 12, height = 10, dpi = 700)

cat("✔ gene_desert_boxplot.pdf saved with human-readable labels\n")
