# clover 0.0.0.9000

* `plot_tRNA_structure()` now centers modification circles, outline circles, and linkage arcs on the visual center of nucleotide characters instead of the text baseline position.

* `plot_tRNA_structure()` overlays modification highlights and circuit linkage arcs on tRNA cloverleaf secondary structure SVGs. Use `structure_organisms()` and `structure_trnas()` to list bundled structures. Base SVGs are generated offline from gtRNAdb data via R2R.

* `structure_to_png()` converts a tRNA structure SVG to PNG format. Requires the rsvg package.

* `plot_tRNA_structure()` linkage arcs now route outward from the structure centroid instead of using a fixed perpendicular offset, avoiding arcs that cut through the interior. Overlapping arcs are assigned to separate lanes for visual clarity. Arc color encodes sign (blue for exclusive, vermillion for co-occurring) and stroke width encodes magnitude of the value. The legend updates to show bidirectional entries when both positive and negative values are present.

* Column naming is now standardized across the package: charging data and DESeq2 results use `ref` instead of `tRNA` for the tRNA reference name, and Sprinzl coordinates use `pos` instead of `seq_index` for the 1-based position in the tRNA body. The default `lab_col` parameter in `plot_volcano()`, `plot_abundance_charging()`, and `tabulate_deseq()` changed from `"tRNA"` to `"ref"`. `read_bcerror()` now returns `ref` as character instead of factor.

* New color palette functions `aa_colors()` and `charging_colors()` provide named color vectors for amino acids and tRNA charging states.

* New statistical utility functions `calc_fold_change()`, `cohens_d()`, `propagate_error_ratio()`, and `propagate_error_diff()` for common tRNA analysis calculations.

* New `trna_regions()` returns a named list mapping canonical tRNA structural region names to Sprinzl position integers.

* `aggregate_or_isodecoder()` collapses per-gene odds ratios to isodecoder level by averaging across gene copies.

* `build_or_network()` constructs a tidygraph network from pairwise odds ratio or ROR data, with node centrality metrics.

* `calculate_rewiring_scores()` summarizes per-isodecoder rewiring magnitude from a ROR matrix.

* `clean_odds_ratios()` prepares odds ratio data for downstream analysis by capping infinite log odds ratio values.

* `compute_ror_isodecoder()` compares isodecoder-level odds ratios between two conditions with z-score significance testing.

* `perform_pcoa()` runs classical multidimensional scaling on a rewiring matrix for dimensionality reduction.

* `plot_arc_diagram()` creates a circular arc diagram from a tidygraph network built by `build_or_network()`.

* `plot_mod_heatmap()` gains new parameters for annotated heatmaps: `label_col` overlays text labels on tiles, `highlight_col` adds dot markers at selected cells, and `group_col` enables group-aware clustering with divider lines between groups. A `caption` parameter adds explanatory text below the plot.

* `plot_mod_landscape()` creates stacked multi-metric profile plots along the tRNA sequence, with optional structural region shading and Sprinzl position secondary axis. Uses patchwork for panel layout.

* `plot_pcoa_rewiring()` creates a scatter plot of PCoA coordinates colored by rewiring magnitude and labeled with top isodecoders.

* `prepare_rewiring_matrix()` builds a wide matrix from isodecoder-level relative odds ratios suitable for PCoA analysis.

* `plot_abundance_charging()` creates a scatter plot comparing tRNA abundance changes (from DESeq2) with charging ratio changes, with significant points colored by quadrant and labeled with ggrepel.

* `plot_bcerror_profile()` plots per-position base-calling error rates as a faceted line plot, with optional modification position overlay.

* `plot_charging_diffs()` creates a dot plot with error bars showing per-tRNA charging ratio differences between conditions.

* `plot_chord_or()` and `plot_chord_ror()` now convert positions from seq_index to Sprinzl labels when `sprinzl_coords` is provided, display all tRNA positions as sectors for structural context, and add annotation rings for structural region, reference nucleotide, and modification positions (via new `mods` parameter). Sectors now have equal widths for consistent visual comparison, a chord color legend is displayed, and default significance cutoffs are tighter (`or_cutoff = 1.0`, `p_cutoff = 0.01`, `min_obs = 100`) to reduce visual clutter.

* `compute_charging_diffs()` compares per-tRNA charging ratios between two conditions, returning mean ratios, standard errors, and the between-condition difference with propagated SE.

* `compute_odds_ratios()` now uses a C++ implementation (via cpp11) for the pairwise Fisher's exact test inner loop, dramatically improving performance on large datasets. The odds ratio is now computed as the sample odds ratio with Haldane correction for zero cells, rather than the conditional MLE.

* `fetch_modomics_mods()` downloads tRNA modification annotations from the MODOMICS database and maps them onto reference sequences using pairwise alignment (#11).

* `modomics_mods()` maps MODOMICS tRNA modifications onto reference sequences using bundled data, eliminating the need for internet access. Use `modomics_organisms()` to list organisms with cached data. Falls back to `fetch_modomics_mods()` for unsupported organisms (#11).

* `plot_volcano()` creates a labeled volcano plot from `tidy_deseq_results()` output, with significant points highlighted and labeled using ggrepel.

* `tabulate_deseq()` creates a formatted gt table of the top significant tRNAs from `tidy_deseq_results()` output, sorted by p-value. The table now includes zebra striping, search, column sorting, and pagination.

* Initial CRAN submission.
