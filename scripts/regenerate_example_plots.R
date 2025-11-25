#!/usr/bin/env Rscript
# Regenerate all example plots in man/figures/
# Run from project root with: Rscript scripts/regenerate_example_plots.R

library(devtools)
library(ggplot2)
library(dplyr)

message("Loading package...")
load_all(".")

message("Loading yeast example data...")
bcerr <- read_bcerror(clover_example("yeast/grande.bcerr.tsv.gz")) |>
  add_global_coords("sacCer") |>
  filter(grepl("^nuc-", ref))

# Plot 1: Heatmap split by tRNA type
message("Generating heatmap-split-example.png...")
p_split <- plot_bcerror_heatmap(bcerr, value="error_rate", show_regions=TRUE)
ggsave("man/figures/heatmap-split-example.png", p_split,
       width=12, height=10, dpi=150, bg="white")

# Plot 2: Unified heatmap
message("Generating heatmap-example.png...")
p_single <- plot_bcerror_heatmap(bcerr, value="error_rate",
                                  show_regions=TRUE, split_by_type=FALSE)
ggsave("man/figures/heatmap-example.png", p_single,
       width=12, height=8, dpi=150, bg="white")

# Plot 3: Single tRNA structure
message("Generating structure-single-example.png...")
p_single_trna <- plot_trna_structure(
  bcerr,
  trna_id = "nuc-tRNA-Ala-AGC-1-1",
  fill = "error_rate",
  modomics = TRUE,
  organism = "sacCer"
)
ggsave("man/figures/structure-single-example.png", p_single_trna,
       width=8, height=8, dpi=150, bg="white")

# Plot 4: Aggregated consensus structure
message("Generating structure-aggregated-example.png...")
p_agg <- plot_trna_structure(bcerr, fill="error_rate")
ggsave("man/figures/structure-aggregated-example.png", p_agg,
       width=8, height=8, dpi=150, bg="white")

# Plot 5: Comparison of grande vs petite
message("Generating structure-comparison-example.png...")
bcerr_grande <- read_bcerror(clover_example("yeast/grande.bcerr.tsv.gz")) |>
  add_global_coords("sacCer") |>
  mutate(condition = "grande")

bcerr_petite <- read_bcerror(clover_example("yeast/petite.bcerr.tsv.gz")) |>
  add_global_coords("sacCer") |>
  mutate(condition = "petite")

bcerr_both <- bind_rows(bcerr_grande, bcerr_petite)

p_comp <- plot_trna_structure(
  bcerr_both,
  compare = c("grande", "petite"),
  fill = "error_rate"
)
ggsave("man/figures/structure-comparison-example.png", p_comp,
       width=12, height=8, dpi=150, bg="white")

message("All plots regenerated successfully!")
