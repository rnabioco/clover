# Package index

## Read

- [`read_bcerror()`](https://rnabioco.github.io/clover/reference/read_bcerror.md)
  : Read base-calling error ("bcerror") TSV files.
- [`read_charging()`](https://rnabioco.github.io/clover/reference/read_charging.md)
  : Read a charging CPM file.
- [`read_charging_multi()`](https://rnabioco.github.io/clover/reference/read_charging_multi.md)
  : Read charging CPM files for multiple samples.
- [`read_counts()`](https://rnabioco.github.io/clover/reference/read_counts.md)
  : Read counts from file
- [`read_fasta()`](https://rnabioco.github.io/clover/reference/read_fasta.md)
  : Read FASTA reference
- [`read_mod_annotations()`](https://rnabioco.github.io/clover/reference/read_mod_annotations.md)
  : Read modifications file
- [`read_odds_ratios()`](https://rnabioco.github.io/clover/reference/read_odds_ratios.md)
  : Read an odds ratios file.
- [`read_odds_ratios_multi()`](https://rnabioco.github.io/clover/reference/read_odds_ratios_multi.md)
  : Read odds ratio files for multiple samples.
- [`read_pipeline_config()`](https://rnabioco.github.io/clover/reference/read_pipeline_config.md)
  : Read a Snakemake pipeline configuration file.
- [`read_pipeline_results()`](https://rnabioco.github.io/clover/reference/read_pipeline_results.md)
  : Read pipeline results for all samples.
- [`read_sprinzl_coords()`](https://rnabioco.github.io/clover/reference/read_sprinzl_coords.md)
  : Read Sprinzl coordinates from a global coordinates TSV file.

## Create

- [`create_clover()`](https://rnabioco.github.io/clover/reference/create_clover.md)
  : Create a SummarizedExperiment from pipeline output.

## Analysis

- [`abundance_count_matrix()`](https://rnabioco.github.io/clover/reference/abundance_count_matrix.md)
  : Build a tRNA abundance count matrix from charging data.
- [`build_coldata()`](https://rnabioco.github.io/clover/reference/build_coldata.md)
  : Build column metadata for DESeq2.
- [`charging_count_matrix()`](https://rnabioco.github.io/clover/reference/charging_count_matrix.md)
  : Build a charging count matrix for differential charging analysis.
- [`compute_charging_diffs()`](https://rnabioco.github.io/clover/reference/compute_charging_diffs.md)
  : Compute charging ratio differences between conditions.
- [`compute_odds_ratios()`](https://rnabioco.github.io/clover/reference/compute_odds_ratios.md)
  : Compute pairwise modification co-occurrence odds ratios.
- [`compute_ror()`](https://rnabioco.github.io/clover/reference/compute_ror.md)
  : Compute ratio of odds ratios between conditions.
- [`run_deseq()`](https://rnabioco.github.io/clover/reference/run_deseq.md)
  : Run DESeq2 differential analysis.
- [`tabulate_deseq()`](https://rnabioco.github.io/clover/reference/tabulate_deseq.md)
  : Tabulate top DESeq2 differential expression results.
- [`tidy_deseq_results()`](https://rnabioco.github.io/clover/reference/tidy_deseq_results.md)
  : Tidy DESeq2 results into a tibble.

## Plot

- [`plot_abundance_charging()`](https://rnabioco.github.io/clover/reference/plot_abundance_charging.md)
  : Plot abundance changes versus charging ratio changes.
- [`plot_bcerror_profile()`](https://rnabioco.github.io/clover/reference/plot_bcerror_profile.md)
  : Plot per-position base-calling error profiles.
- [`plot_charging_diffs()`](https://rnabioco.github.io/clover/reference/plot_charging_diffs.md)
  : Plot per-tRNA charging ratio differences.
- [`plot_chord_or()`](https://rnabioco.github.io/clover/reference/plot_chord_or.md)
  : Plot a chord diagram of modification co-occurrence.
- [`plot_chord_ror()`](https://rnabioco.github.io/clover/reference/plot_chord_ror.md)
  : Plot a chord diagram of modification rewiring between conditions.
- [`plot_mod_heatmap()`](https://rnabioco.github.io/clover/reference/plot_mod_heatmap.md)
  : Plot a delta-signal modification heatmap.
- [`plot_volcano()`](https://rnabioco.github.io/clover/reference/plot_volcano.md)
  : Plot a volcano plot of differential expression results.

## MODOMICS

- [`modomics_mods()`](https://rnabioco.github.io/clover/reference/modomics_mods.md)
  : Map cached MODOMICS modifications onto reference sequences
- [`modomics_organisms()`](https://rnabioco.github.io/clover/reference/modomics_organisms.md)
  : List organisms with cached MODOMICS data
- [`fetch_modomics_mods()`](https://rnabioco.github.io/clover/reference/fetch_modomics_mods.md)
  : Fetch tRNA modification annotations from MODOMICS

## Utilities

- [`clover_example()`](https://rnabioco.github.io/clover/reference/clover_example.md)
  : Provide working directory for clover example files.
- [`list_pipeline_files()`](https://rnabioco.github.io/clover/reference/list_pipeline_files.md)
  : List pipeline output files for each sample.
- [`order_sprinzl_positions()`](https://rnabioco.github.io/clover/reference/order_sprinzl_positions.md)
  : Order Sprinzl position labels.
