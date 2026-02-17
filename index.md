# clover

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

## Example

clover reads output from the
[aa-tRNA-seq-pipeline](https://github.com/rnabioco/aa-tRNA-seq-pipeline)
and stores the results in a `SummarizedExperiment`.

``` r
library(clover)

# Load pipeline results into a SummarizedExperiment
se <- create_clover(
  config_path = clover_example("ecoli/config.yaml"),
  sample_info = data.frame(
    sample_id = c(
      "wt-15-ctl-01", "wt-15-ctl-02", "wt-15-ctl-03",
      "wt-15-inf-01", "wt-15-inf-02", "wt-15-inf-03"
    ),
    condition = rep(c("control", "infected"), each = 3)
  )
)

# Differential tRNA abundance with DESeq2
dds <- run_deseq(se, design = ~condition)
res <- tidy_deseq_results(dds, contrast = c("condition", "infected", "control"))
plot_volcano(res)

# Base-calling error heatmap with Sprinzl coordinates
sprinzl <- read_sprinzl_coords(
  clover_example("sprinzl/ecoliK12_global_coords.tsv.gz")
)

plot_mod_heatmap(
  S4Vectors::metadata(se)$bcerror,
  sprinzl_coords = sprinzl
)
```

See
[`vignette("clover")`](https://rnabioco.github.io/clover/articles/clover.md)
for a complete walkthrough.

## Related work

- [aa-tRNA-seq-pipeline](https://github.com/rnabioco/aa-tRNA-seq-pipeline)
  is the Snakemake pipeline that generates the data clover analyzes.
- [R2easyR](https://github.com/JPSieg/R2easyR) visualizes structure
  probing signals on RNA secondary structure diagrams.
- [nanoblot](https://github.com/SamDeMario-lab/NanoBlot) facilitates
  visualization of nanopore sequencing data, including a “virtual gel”
  plot.
