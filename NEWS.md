# clover 0.0.0.9000

* `plot_abundance_charging()` creates a scatter plot comparing tRNA abundance changes (from DESeq2) with charging ratio changes, with significant points colored by quadrant and labeled with ggrepel.

* `plot_bcerror_profile()` plots per-position base-calling error rates as a faceted line plot, with optional modification position overlay.

* `plot_charging_diffs()` creates a dot plot with error bars showing per-tRNA charging ratio differences between conditions.

* `plot_chord_or()` and `plot_chord_ror()` now convert positions from seq_index to Sprinzl labels when `sprinzl_coords` is provided, display all tRNA positions as sectors for structural context, and add annotation rings for structural region, reference nucleotide, and modification positions (via new `mods` parameter). Sectors now have equal widths for consistent visual comparison, a chord color legend is displayed, and default significance cutoffs are tighter (`or_cutoff = 1.0`, `p_cutoff = 0.01`, `min_obs = 100`) to reduce visual clutter.

* `compute_charging_diffs()` compares per-tRNA charging ratios between two conditions, returning mean ratios, standard errors, and the between-condition difference with propagated SE.

* `compute_odds_ratios()` now uses a C++ implementation (via cpp11) for the pairwise Fisher's exact test inner loop, dramatically improving performance on large datasets. The odds ratio is now computed as the sample odds ratio with Haldane correction for zero cells, rather than the conditional MLE.

* `fetch_modomics_mods()` downloads tRNA modification annotations from the MODOMICS database and maps them onto reference sequences using pairwise alignment (#11).

* `modomics_mods()` maps MODOMICS tRNA modifications onto reference sequences using bundled data, eliminating the need for internet access. Use `modomics_organisms()` to list organisms with cached data. Falls back to `fetch_modomics_mods()` for unsupported organisms (#11).

* `plot_volcano()` creates a labeled volcano plot from `tidy_deseq_results()` output, with significant points highlighted and labeled using ggrepel.

* Initial CRAN submission.
