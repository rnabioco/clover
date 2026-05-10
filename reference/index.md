# Package index

## Read

- [`read_bcerror()`](https://rnabioco.github.io/clover/reference/read_bcerror.md)
  : Read base-calling error ("bcerror") TSV files.
- [`read_bedmethyl()`](https://rnabioco.github.io/clover/reference/read_bedmethyl.md)
  : Read a bedMethyl file from modkit pileup
- [`read_charging()`](https://rnabioco.github.io/clover/reference/read_charging.md)
  : Read a charging CPM file.
- [`read_charging_multi()`](https://rnabioco.github.io/clover/reference/read_charging_multi.md)
  : Read charging CPM files for multiple samples.
- [`read_counts()`](https://rnabioco.github.io/clover/reference/read_counts.md)
  : Read a counts TSV file.
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
- [`aggregate_or_isodecoder()`](https://rnabioco.github.io/clover/reference/aggregate_or_isodecoder.md)
  : Aggregate odds ratios to isodecoder level.
- [`build_coldata()`](https://rnabioco.github.io/clover/reference/build_coldata.md)
  : Build column metadata for DESeq2.
- [`calculate_rewiring_scores()`](https://rnabioco.github.io/clover/reference/calculate_rewiring_scores.md)
  : Calculate rewiring scores from a ROR matrix.
- [`charging_count_matrix()`](https://rnabioco.github.io/clover/reference/charging_count_matrix.md)
  : Build a charging count matrix for differential charging analysis.
- [`clean_odds_ratios()`](https://rnabioco.github.io/clover/reference/clean_odds_ratios.md)
  : Clean odds ratio data.
- [`compute_bcerror_delta()`](https://rnabioco.github.io/clover/reference/compute_bcerror_delta.md)
  : Compute per-position delta between two conditions.
- [`compute_charging_diffs()`](https://rnabioco.github.io/clover/reference/compute_charging_diffs.md)
  : Compute charging ratio differences between conditions.
- [`filter_linkages()`](https://rnabioco.github.io/clover/reference/filter_linkages.md)
  : Filter odds ratios for structure linkage arcs.
- [`compute_odds_ratios()`](https://rnabioco.github.io/clover/reference/compute_odds_ratios.md)
  : Compute pairwise modification co-occurrence odds ratios.
- [`summarize_mod_calls()`](https://rnabioco.github.io/clover/reference/summarize_mod_calls.md)
  : Summarize per-position modification frequency from modkit extract
- [`compute_ror()`](https://rnabioco.github.io/clover/reference/compute_ror.md)
  : Compute ratio of odds ratios between conditions.
- [`compute_ror_isodecoder()`](https://rnabioco.github.io/clover/reference/compute_ror_isodecoder.md)
  : Compute relative odds ratio at isodecoder level.
- [`perform_pcoa()`](https://rnabioco.github.io/clover/reference/perform_pcoa.md)
  : Perform PCoA on a rewiring matrix.
- [`prepare_rewiring_matrix()`](https://rnabioco.github.io/clover/reference/prepare_rewiring_matrix.md)
  : Prepare a matrix for rewiring PCoA analysis.
- [`run_deseq()`](https://rnabioco.github.io/clover/reference/run_deseq.md)
  : Run DESeq2 differential analysis.
- [`tabulate_deseq()`](https://rnabioco.github.io/clover/reference/tabulate_deseq.md)
  : Tabulate top DESeq2 differential expression results.
- [`tidy_deseq_results()`](https://rnabioco.github.io/clover/reference/tidy_deseq_results.md)
  : Tidy DESeq2 results into a tibble.

## Plot

- [`prep_mod_heatmap()`](https://rnabioco.github.io/clover/reference/prep_mod_heatmap.md)
  :

  Prepare data for
  [`plot_mod_heatmap()`](https://rnabioco.github.io/clover/reference/plot_mod_heatmap.md).

- [`plot_abundance_charging()`](https://rnabioco.github.io/clover/reference/plot_abundance_charging.md)
  : Plot abundance changes versus charging ratio changes.

- [`plot_bcerror_profile()`](https://rnabioco.github.io/clover/reference/plot_bcerror_profile.md)
  : Plot per-position base-calling error profiles.

- [`plot_charging_diffs()`](https://rnabioco.github.io/clover/reference/plot_charging_diffs.md)
  : Plot per-tRNA charging ratio differences.

- [`plot_charging_ratios()`](https://rnabioco.github.io/clover/reference/plot_charging_ratios.md)
  : Plot per-tRNA charging ratios.

- [`plot_chord_or()`](https://rnabioco.github.io/clover/reference/plot_chord_or.md)
  : Plot a chord diagram of modification co-occurrence.

- [`plot_chord_ror()`](https://rnabioco.github.io/clover/reference/plot_chord_ror.md)
  : Plot a chord diagram of modification rewiring between conditions.

- [`plot_mod_heatmap()`](https://rnabioco.github.io/clover/reference/plot_mod_heatmap.md)
  : Plot a delta-signal modification heatmap.

- [`plot_mod_landscape()`](https://rnabioco.github.io/clover/reference/plot_mod_landscape.md)
  : Plot per-tRNA modification landscape profiles.

- [`plot_pcoa_rewiring()`](https://rnabioco.github.io/clover/reference/plot_pcoa_rewiring.md)
  : Plot PCoA of tRNA rewiring scores.

- [`plot_volcano()`](https://rnabioco.github.io/clover/reference/plot_volcano.md)
  : Plot a volcano plot of differential expression results.

## Network

- [`build_or_network()`](https://rnabioco.github.io/clover/reference/build_or_network.md)
  : Build a network from odds ratio data.
- [`plot_arc_diagram()`](https://rnabioco.github.io/clover/reference/plot_arc_diagram.md)
  : Plot a network as an arc diagram.

## Structure

- [`default_mod_palette()`](https://rnabioco.github.io/clover/reference/default_mod_palette.md)
  : Default modification color palette
- [`plot_tRNA_structure()`](https://rnabioco.github.io/clover/reference/plot_tRNA_structure.md)
  : Plot tRNA secondary structure with modifications and linkages
- [`structure_html()`](https://rnabioco.github.io/clover/reference/structure_html.md)
  : Embed a tRNA structure SVG as centered HTML
- [`structure_organisms()`](https://rnabioco.github.io/clover/reference/structure_organisms.md)
  : List organisms with bundled tRNA structure SVGs
- [`structure_to_png()`](https://rnabioco.github.io/clover/reference/structure_to_png.md)
  : Convert a tRNA structure SVG to PNG
- [`structure_trnas()`](https://rnabioco.github.io/clover/reference/structure_trnas.md)
  : List available tRNA structures for an organism

## MODOMICS

- [`modomics_mods()`](https://rnabioco.github.io/clover/reference/modomics_mods.md)
  : Map cached MODOMICS modifications onto reference sequences
- [`modomics_organisms()`](https://rnabioco.github.io/clover/reference/modomics_organisms.md)
  : List organisms with cached MODOMICS data
- [`fetch_modomics_mods()`](https://rnabioco.github.io/clover/reference/fetch_modomics_mods.md)
  : Fetch tRNA modification annotations from MODOMICS

## Identity elements

- [`identity_elements()`](https://rnabioco.github.io/clover/reference/identity_elements.md)
  : Retrieve tRNA aminoacylation identity elements
- [`identity_organisms()`](https://rnabioco.github.io/clover/reference/identity_organisms.md)
  : List supported organisms for identity elements
- [`map_identity_to_trna()`](https://rnabioco.github.io/clover/reference/map_identity_to_trna.md)
  : Map identity elements to tRNA sequence positions
- [`plot_identity_panel()`](https://rnabioco.github.io/clover/reference/plot_identity_panel.md)
  : Plot identity elements for multiple tRNAs in a grid
- [`plot_identity_structure()`](https://rnabioco.github.io/clover/reference/plot_identity_structure.md)
  : Plot tRNA structure with identity element overlays
- [`tertiary_contacts()`](https://rnabioco.github.io/clover/reference/tertiary_contacts.md)
  : Retrieve canonical tRNA tertiary contacts

## Statistics

- [`calc_fold_change()`](https://rnabioco.github.io/clover/reference/calc_fold_change.md)
  : Calculate log2 fold change with pseudocount.
- [`cohens_d()`](https://rnabioco.github.io/clover/reference/cohens_d.md)
  : Calculate Cohen's d effect size.
- [`propagate_error_diff()`](https://rnabioco.github.io/clover/reference/propagate_error_diff.md)
  : Propagate standard error for a difference.
- [`propagate_error_ratio()`](https://rnabioco.github.io/clover/reference/propagate_error_ratio.md)
  : Propagate standard error for a ratio.

## Utilities

- [`aa_colors()`](https://rnabioco.github.io/clover/reference/aa_colors.md)
  : Amino acid color palette.
- [`charging_colors()`](https://rnabioco.github.io/clover/reference/charging_colors.md)
  : Charging status color palette.
- [`clover_example()`](https://rnabioco.github.io/clover/reference/clover_example.md)
  : Provide working directory for clover example files.
- [`dna_to_rna_anticodon()`](https://rnabioco.github.io/clover/reference/dna_to_rna_anticodon.md)
  : Convert DNA anticodon to RNA in tRNA name strings
- [`shorten_trna_names()`](https://rnabioco.github.io/clover/reference/shorten_trna_names.md)
  : Shorten tRNA name strings for display.
- [`list_pipeline_files()`](https://rnabioco.github.io/clover/reference/list_pipeline_files.md)
  : List pipeline output files for each sample.
- [`order_sprinzl_positions()`](https://rnabioco.github.io/clover/reference/order_sprinzl_positions.md)
  : Order Sprinzl position labels.
- [`trna_regions()`](https://rnabioco.github.io/clover/reference/trna_regions.md)
  : tRNA structural regions.
