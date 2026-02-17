# Project Tracker

## 2026-02-16: E. coli test data, reader updates, create_clover()

### Completed

- Added E. coli wt-15 test data (6 samples: 3 ctl + 3 inf) to `inst/extdata/ecoli/`
  - charging, bcerror (5 tRNAs subset), align_stats, odds_ratios, reference FASTA
  - Pipeline directory structure: `summary/tables/{sample_id}/`
  - config.yaml + headerless samples.tsv matching real pipeline format
- Rewrote `read_bcerror()` for new pipeline column format (Reference, Position, Spanning_Reads, etc.)
- Updated `read_pipeline_config()` to support `output_directory` key and headerless samples.tsv
- Made `read_pipeline_results()` skip missing files gracefully
- Made `.resolve_path()` robust to non-character YAML values
- Added `compute_odds_ratios()` function in `R/read-charging.R`
- Implemented `create_clover()` in `R/clover-se.R` replacing skeleton `CloverSE()`
  - Builds SummarizedExperiment with counts assay, colData, rowData, metadata
  - Loads bcerror and odds_ratios into metadata
  - Reads reference FASTA and adds sequence lengths to rowData
- Added S4Vectors to DESCRIPTION imports
- Updated globals.R with new variable names
- Updated all tests (136 pass, 0 fail)
- Removed old yeast bcerror files (grande.bcerr.tsv.gz, petite.bcerr.tsv.gz)

## 2026-02-16: Vignette and chord diagram fixes

### Completed

- Created `vignettes/ecoli-phage.qmd` demonstrating full analysis workflow:
  - `create_clover()` data loading
  - Differential abundance with DESeq2 (volcano plot)
  - Charging ratio analysis (dot plot)
  - Base-calling error profiles (line plots)
  - Error rate difference heatmap with Sprinzl coordinates
  - Modification co-occurrence chord diagrams (`plot_chord_or()`)
  - Modification rewiring chord diagrams (`compute_ror()` + `plot_chord_ror()`)
  - Step-by-step data loading alternative
- Fixed chord diagrams rendering blank: root cause was `Inf` values in
  `log_odds_ratio` from Fisher's exact test (when one cell in 2x2 table is 0).
  Added `is.finite()` filter in `plot_chord_or()`, `plot_chord_ror()`, and
  `compute_ror()`.
- Added `.map_positions()` helper to convert numeric seq_index positions to
  Sprinzl labels when `sprinzl_coords` are provided
- Fixed `.region_colors()` palette to match actual region names in sprinzl data
  (e.g., "anticodon-stem" not "AC-stem", "variable-region" not "variable-loop")
- Refactored `.draw_chord()` as shared helper for consistent chord diagram styling
- All 136 tests pass

### Known issues

- `R CMD build` / `devtools::check()` fails due to pre-existing pixi GCC issue
  (file with `TOOLS=` in name confuses `cp`). Not related to code changes.
  `devtools::load_all()` + `devtools::test()` work fine.

### Next steps

- Fix `R CMD build` issue (possibly update `.Rbuildignore` or pixi env)
- Add integration tests for DESeq2 with real E. coli data
- Add plot tests with real E. coli data
- Update README.Rmd to use new E. coli data instead of yeast
