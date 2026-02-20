# Changelog

## clover 0.0.0.9000

- [`plot_tRNA_structure()`](https://rnabioco.github.io/clover/reference/plot_tRNA_structure.md)
  now centers modification circles, outline circles, and linkage arcs on
  the visual center of nucleotide characters instead of the text
  baseline position.

- [`plot_tRNA_structure()`](https://rnabioco.github.io/clover/reference/plot_tRNA_structure.md)
  overlays modification highlights and circuit linkage arcs on tRNA
  cloverleaf secondary structure SVGs. Use
  [`structure_organisms()`](https://rnabioco.github.io/clover/reference/structure_organisms.md)
  and
  [`structure_trnas()`](https://rnabioco.github.io/clover/reference/structure_trnas.md)
  to list bundled structures. Base SVGs are generated offline from
  gtRNAdb data via R2R.

- [`structure_to_png()`](https://rnabioco.github.io/clover/reference/structure_to_png.md)
  converts a tRNA structure SVG to PNG format. Requires the rsvg
  package.

- [`plot_tRNA_structure()`](https://rnabioco.github.io/clover/reference/plot_tRNA_structure.md)
  linkage arcs now route outward from the structure centroid instead of
  using a fixed perpendicular offset, avoiding arcs that cut through the
  interior. Overlapping arcs are assigned to separate lanes for visual
  clarity. Arc color encodes sign (blue for exclusive, vermillion for
  co-occurring) and stroke width encodes magnitude of the value. The
  legend updates to show bidirectional entries when both positive and
  negative values are present.

- Column naming is now standardized across the package: charging data
  and DESeq2 results use `ref` instead of `tRNA` for the tRNA reference
  name, and Sprinzl coordinates use `pos` instead of `seq_index` for the
  1-based position in the tRNA body. The default `lab_col` parameter in
  [`plot_volcano()`](https://rnabioco.github.io/clover/reference/plot_volcano.md),
  [`plot_abundance_charging()`](https://rnabioco.github.io/clover/reference/plot_abundance_charging.md),
  and
  [`tabulate_deseq()`](https://rnabioco.github.io/clover/reference/tabulate_deseq.md)
  changed from `"tRNA"` to `"ref"`.
  [`read_bcerror()`](https://rnabioco.github.io/clover/reference/read_bcerror.md)
  now returns `ref` as character instead of factor.

- New color palette functions
  [`aa_colors()`](https://rnabioco.github.io/clover/reference/aa_colors.md)
  and
  [`charging_colors()`](https://rnabioco.github.io/clover/reference/charging_colors.md)
  provide named color vectors for amino acids and tRNA charging states.

- New statistical utility functions
  [`calc_fold_change()`](https://rnabioco.github.io/clover/reference/calc_fold_change.md),
  [`cohens_d()`](https://rnabioco.github.io/clover/reference/cohens_d.md),
  [`propagate_error_ratio()`](https://rnabioco.github.io/clover/reference/propagate_error_ratio.md),
  and
  [`propagate_error_diff()`](https://rnabioco.github.io/clover/reference/propagate_error_diff.md)
  for common tRNA analysis calculations.

- New
  [`trna_regions()`](https://rnabioco.github.io/clover/reference/trna_regions.md)
  returns a named list mapping canonical tRNA structural region names to
  Sprinzl position integers.

- [`aggregate_or_isodecoder()`](https://rnabioco.github.io/clover/reference/aggregate_or_isodecoder.md)
  collapses per-gene odds ratios to isodecoder level by averaging across
  gene copies.

- [`build_or_network()`](https://rnabioco.github.io/clover/reference/build_or_network.md)
  constructs a tidygraph network from pairwise odds ratio or ROR data,
  with node centrality metrics.

- [`calculate_rewiring_scores()`](https://rnabioco.github.io/clover/reference/calculate_rewiring_scores.md)
  summarizes per-isodecoder rewiring magnitude from a ROR matrix.

- [`clean_odds_ratios()`](https://rnabioco.github.io/clover/reference/clean_odds_ratios.md)
  prepares odds ratio data for downstream analysis by capping infinite
  log odds ratio values.

- [`compute_ror_isodecoder()`](https://rnabioco.github.io/clover/reference/compute_ror_isodecoder.md)
  compares isodecoder-level odds ratios between two conditions with
  z-score significance testing.

- [`perform_pcoa()`](https://rnabioco.github.io/clover/reference/perform_pcoa.md)
  runs classical multidimensional scaling on a rewiring matrix for
  dimensionality reduction.

- [`plot_arc_diagram()`](https://rnabioco.github.io/clover/reference/plot_arc_diagram.md)
  creates a circular arc diagram from a tidygraph network built by
  [`build_or_network()`](https://rnabioco.github.io/clover/reference/build_or_network.md).

- [`plot_mod_heatmap()`](https://rnabioco.github.io/clover/reference/plot_mod_heatmap.md)
  gains new parameters for annotated heatmaps: `label_col` overlays text
  labels on tiles, `highlight_col` adds dot markers at selected cells,
  and `group_col` enables group-aware clustering with divider lines
  between groups. A `caption` parameter adds explanatory text below the
  plot.

- [`plot_mod_landscape()`](https://rnabioco.github.io/clover/reference/plot_mod_landscape.md)
  creates stacked multi-metric profile plots along the tRNA sequence,
  with optional structural region shading and Sprinzl position secondary
  axis. Uses patchwork for panel layout.

- [`plot_pcoa_rewiring()`](https://rnabioco.github.io/clover/reference/plot_pcoa_rewiring.md)
  creates a scatter plot of PCoA coordinates colored by rewiring
  magnitude and labeled with top isodecoders.

- [`prepare_rewiring_matrix()`](https://rnabioco.github.io/clover/reference/prepare_rewiring_matrix.md)
  builds a wide matrix from isodecoder-level relative odds ratios
  suitable for PCoA analysis.

- [`plot_abundance_charging()`](https://rnabioco.github.io/clover/reference/plot_abundance_charging.md)
  creates a scatter plot comparing tRNA abundance changes (from DESeq2)
  with charging ratio changes, with significant points colored by
  quadrant and labeled with ggrepel.

- [`plot_bcerror_profile()`](https://rnabioco.github.io/clover/reference/plot_bcerror_profile.md)
  plots per-position base-calling error rates as a faceted line plot,
  with optional modification position overlay.

- [`plot_charging_diffs()`](https://rnabioco.github.io/clover/reference/plot_charging_diffs.md)
  creates a dot plot with error bars showing per-tRNA charging ratio
  differences between conditions.

- [`plot_chord_or()`](https://rnabioco.github.io/clover/reference/plot_chord_or.md)
  and
  [`plot_chord_ror()`](https://rnabioco.github.io/clover/reference/plot_chord_ror.md)
  now convert positions from seq_index to Sprinzl labels when
  `sprinzl_coords` is provided, display all tRNA positions as sectors
  for structural context, and add annotation rings for structural
  region, reference nucleotide, and modification positions (via new
  `mods` parameter). Sectors now have equal widths for consistent visual
  comparison, a chord color legend is displayed, and default
  significance cutoffs are tighter (`or_cutoff = 1.0`,
  `p_cutoff = 0.01`, `min_obs = 100`) to reduce visual clutter.

- [`compute_charging_diffs()`](https://rnabioco.github.io/clover/reference/compute_charging_diffs.md)
  compares per-tRNA charging ratios between two conditions, returning
  mean ratios, standard errors, and the between-condition difference
  with propagated SE.

- [`compute_odds_ratios()`](https://rnabioco.github.io/clover/reference/compute_odds_ratios.md)
  now uses a C++ implementation (via cpp11) for the pairwise Fisher’s
  exact test inner loop, dramatically improving performance on large
  datasets. The odds ratio is now computed as the sample odds ratio with
  Haldane correction for zero cells, rather than the conditional MLE.

- [`fetch_modomics_mods()`](https://rnabioco.github.io/clover/reference/fetch_modomics_mods.md)
  downloads tRNA modification annotations from the MODOMICS database and
  maps them onto reference sequences using pairwise alignment
  ([\#11](https://github.com/rnabioco/clover/issues/11)).

- [`modomics_mods()`](https://rnabioco.github.io/clover/reference/modomics_mods.md)
  maps MODOMICS tRNA modifications onto reference sequences using
  bundled data, eliminating the need for internet access. Use
  [`modomics_organisms()`](https://rnabioco.github.io/clover/reference/modomics_organisms.md)
  to list organisms with cached data. Falls back to
  [`fetch_modomics_mods()`](https://rnabioco.github.io/clover/reference/fetch_modomics_mods.md)
  for unsupported organisms
  ([\#11](https://github.com/rnabioco/clover/issues/11)).

- [`plot_volcano()`](https://rnabioco.github.io/clover/reference/plot_volcano.md)
  creates a labeled volcano plot from
  [`tidy_deseq_results()`](https://rnabioco.github.io/clover/reference/tidy_deseq_results.md)
  output, with significant points highlighted and labeled using ggrepel.

- [`tabulate_deseq()`](https://rnabioco.github.io/clover/reference/tabulate_deseq.md)
  creates a formatted gt table of the top significant tRNAs from
  [`tidy_deseq_results()`](https://rnabioco.github.io/clover/reference/tidy_deseq_results.md)
  output, sorted by p-value. The table now includes zebra striping,
  search, column sorting, and pagination.

- Initial CRAN submission.
