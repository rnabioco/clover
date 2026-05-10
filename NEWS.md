# clover 0.0.0.9000

* `read_bedmethyl()` reads per-position modification percentages from bedMethyl files produced by `modkit pileup`, enabling integration of nanopore modification data with the heatmap workflow.

* `summarize_mod_calls()` summarizes per-read modification calls from `modkit extract` output into per-position modification frequencies, ready for delta computation and heatmap visualization.

* `plot_abundance_charging()` gains `shorten`, `source_col`, and `error_bars` parameters. Labels are auto-shortened via `shorten_trna_names()`, an optional faceting column separates host and phage tRNAs, and horizontal error bars show `lfcSE` on significant points when present.

* `plot_charging_diffs()` gains `source_col`, `label_col`, and `shorten` parameters. Labels are auto-shortened via `shorten_trna_names()` and an optional faceting column separates host and phage tRNAs.

* `plot_charging_ratios()` creates a box-and-jitter plot of raw charging ratios grouped by condition, with optional faceting and auto-shortened tRNA labels.

* `shorten_trna_names()` extracts the tRNA label-shortening logic into a reusable utility function, stripping source prefixes, the `tRNA-` prefix, and gene-copy suffixes.

* `compute_bcerror_delta()` computes per-position differences in base-calling error rates between two conditions from a summarized bcerror tibble.

* `prep_mod_heatmap()` prepares bcerror delta data for `plot_mod_heatmap()` by joining Sprinzl coordinates, annotating known modifications, and shortening tRNA labels.

* `identity_elements()` returns experimentally validated tRNA aminoacylation identity elements (determinants and antideterminants) for a given organism, based on Giege & Eriani (2023). Use `identity_organisms()` to list supported organisms.

* `map_identity_to_trna()` converts Sprinzl-numbered identity elements to 1-based sequence positions for a specific tRNA, enabling overlay on `plot_tRNA_structure()` via the `outlines` parameter.

* `plot_identity_panel()` arranges multiple tRNA cloverleaf structures side by side in a single SVG, each annotated with aminoacylation identity elements colored by strength (red = strong, blue = weak).

* `plot_identity_structure()` generates a single tRNA cloverleaf structure SVG with identity element outlines, automatically mapping Sprinzl positions to sequence coordinates.

* `structure_html()` wraps a tRNA structure SVG in a centering `<div>` and returns an `htmltools::HTML` object, simplifying embedding in R Markdown and Quarto documents.

* `filter_linkages()` filters odds ratio data by p-value, observation count, and log odds ratio magnitude, returning a tibble ready for the `linkages` parameter of `plot_tRNA_structure()`.

* `plot_tRNA_structure()` gains a `layout` argument with `"cloverleaf"` (default) and `"elbow"` forms; the elbow form stacks the acceptor stem and T-arm coaxially at the top with the D-arm extending horizontally and the anticodon arm vertically. Currently bundled for E. coli standard 3-arm tRNAs. `structure_trnas()` accepts the same `layout` argument.

* `plot_tRNA_structure()` now accepts odds ratio tibbles directly as `linkages` input: if a `log_odds_ratio` column is present and `value` is not, it is automatically used as the arc value.

* `plot_tRNA_structure()` now draws a 3' amino acid label (e.g., "Glu") connected by a line to the terminal nucleotide, and position markers every 10 nucleotides around the cloverleaf. Position markers can be disabled with `position_markers = FALSE`.

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

* `plot_mod_heatmap()` gains `cluster_threshold` to filter low-magnitude noise before clustering, `fill_name` and `fill_breaks` to customize the legend title and breaks, and improved caption styling with `plot.caption.position = "plot"`.

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
