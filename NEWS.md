# clover 0.0.0.9000

* `compute_charging_diffs()` compares per-tRNA charging ratios between two conditions, returning mean ratios, standard errors, and the between-condition difference with propagated SE.

* `fetch_modomics_mods()` downloads tRNA modification annotations from the MODOMICS database and maps them onto reference sequences using pairwise alignment (#11).

* `plot_volcano()` creates a labeled volcano plot from `tidy_deseq_results()` output, with significant points highlighted and labeled using ggrepel.

* Initial CRAN submission.
