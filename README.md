
# clover

<!-- badges: start -->

[![R-CMD-check](https://github.com/rnabioco/clover/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/rnabioco/clover/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

clover facilitates analysis and visualization of nanopore tRNA
sequencing data, including differential expression, base-calling error
analysis, and modification co-occurrence networks.

**clover is under active development.** *Caveat emptor*.

## Installation

You can install the development version of clover from
[GitHub](https://github.com/rnabioco/clover) with:

``` r
# install.packages("pak")
pak::pak("rnabioco/clover")
```

## Usage

clover reads output from the
[aa-tRNA-seq-pipeline](https://github.com/rnabioco/aa-tRNA-seq-pipeline)
and stores the results in a `SummarizedExperiment`. The main entry point
is `create_clover()`, which reads a pipeline `config.yaml` and loads
counts, base-calling error rates, and modification co-occurrence data.

``` r
library(clover)

sample_info <- data.frame(
  sample_id = c(
    "wt-15-ctl-01",
    "wt-15-ctl-02",
    "wt-15-ctl-03",
    "wt-15-inf-01",
    "wt-15-inf-02",
    "wt-15-inf-03"
  ),
  condition = rep(c("control", "infected"), each = 3)
)

se <- create_clover(
  config_path = clover_example("ecoli/config.yaml"),
  sample_info = sample_info
)

se
```

## Analysis

### Differential tRNA abundance

clover wraps DESeq2 to test for differential tRNA abundance between
conditions.

``` r
dds <- run_deseq(se, design = ~condition)
res <- tidy_deseq_results(dds, contrast = c("condition", "infected", "control"))
```

### Base-calling error profiles

Error rates per tRNA position reveal modification signatures. clover
provides functions for plotting error profiles and heatmaps annotated
with Sprinzl structural coordinates.

``` r
# Error rate line plot for a specific tRNA
plot_bcerror(
  S4Vectors::metadata(se)$bcerror,
  ref = "host-tRNA-Glu-TTC-1-1"
)

# Heatmap of error rate differences with Sprinzl coordinates
sprinzl <- read_sprinzl_coords(
  clover_example("sprinzl/ecoliK12_global_coords.tsv.gz")
)

plot_mod_heatmap(
  S4Vectors::metadata(se)$bcerror,
  sprinzl_coords = sprinzl
)
```

### Modification co-occurrence

Chord diagrams display pairwise modification co-occurrence (odds ratios)
within a single sample or changes between conditions (ratio of odds
ratios).

``` r
or_data <- S4Vectors::metadata(se)$odds_ratios

# Single-sample chord diagram
or_single <- or_data |>
  dplyr::filter(
    ref == "host-tRNA-Glu-TTC-1-1",
    sample_id == "wt-15-ctl-01"
  )

plot_chord_or(or_single, sprinzl_coords = sprinzl)

# Rewiring between conditions
or_data$condition <- ifelse(
  grepl("ctl", or_data$sample_id),
  "control",
  "infected"
)

ror <- compute_ror(
  or_data,
  numerator = "infected",
  denominator = "control"
)

plot_chord_ror(ror, sprinzl_coords = sprinzl)
```

## Related work

- [aa-tRNA-seq-pipeline](https://github.com/rnabioco/aa-tRNA-seq-pipeline)
  is the Snakemake pipeline that generates the data clover analyzes.
- [R2easyR](https://github.com/JPSieg/R2easyR) visualizes structure
  probing signals on RNA secondary structure diagrams.
- [nanoblot](https://github.com/SamDeMario-lab/NanoBlot) facilitates
  visualization of nanopore sequencing data, including a “virtual gel”
  plot.
