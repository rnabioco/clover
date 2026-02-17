# Changelog

## clover 0.0.0.9000

- [`plot_chord_or()`](https://rnabioco.github.io/clover/reference/plot_chord_or.md)
  and
  [`plot_chord_ror()`](https://rnabioco.github.io/clover/reference/plot_chord_ror.md)
  now convert positions from seq_index to Sprinzl labels when
  `sprinzl_coords` is provided, display all tRNA positions as sectors
  for structural context, and add annotation rings for structural
  region, reference nucleotide, and modification positions (via new
  `mods` parameter).

- [`compute_charging_diffs()`](https://rnabioco.github.io/clover/reference/compute_charging_diffs.md)
  compares per-tRNA charging ratios between two conditions, returning
  mean ratios, standard errors, and the between-condition difference
  with propagated SE.

- [`fetch_modomics_mods()`](https://rnabioco.github.io/clover/reference/fetch_modomics_mods.md)
  downloads tRNA modification annotations from the MODOMICS database and
  maps them onto reference sequences using pairwise alignment
  ([\#11](https://github.com/rnabioco/clover/issues/11)).

- [`plot_volcano()`](https://rnabioco.github.io/clover/reference/plot_volcano.md)
  creates a labeled volcano plot from
  [`tidy_deseq_results()`](https://rnabioco.github.io/clover/reference/tidy_deseq_results.md)
  output, with significant points highlighted and labeled using ggrepel.

- Initial CRAN submission.
