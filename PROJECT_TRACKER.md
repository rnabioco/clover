# Project Tracker

Track development progress and session notes for clover.

## Current Status

**Version**: 0.0.0.9000 (active development)
**Branch**: See `git branch` for current working branch

## Recent Sessions

### 2025-11-25 - Offset×Type Coordinate System Integration

**Work completed:**
- Integrated new offset×type coordinate grouping from tRNAs-in-space
- Replaced single unified coordinate files with grouped files per organism
- Each organism now has 3-7 coordinate files (grouped by offset and type):
  - E. coli: 5 groups (82 tRNAs)
  - Yeast: 4 groups (268 tRNAs)
  - Human: 7 groups (421 tRNAs)
- Added auto-classification: tRNAs automatically matched to correct coordinate group
- Created lookup tables (`{organism}_trna_groups.tsv`) for fast group assignment
- Updated Modomics files with enhanced tRNAs-in-space mappings
- Position 55 now aligns correctly (single global_index per group)
- All 112 tests pass with 0 errors, 0 warnings

**Key API changes:**
- `load_global_coords()` now returns `offset` and `type` columns
- `add_global_coords()` auto-classifies tRNAs (warns about unmatched SeC/mito/iMet)
- `available_coord_groups()` new function to list groups per organism
- `get_global_labels()` returns list when multiple groups present

**Files added/modified:**
- `inst/extdata/coords/` - 16 new grouped coordinate files
- `inst/extdata/coords/*_trna_groups.tsv` - lookup tables (3 files)
- `inst/extdata/modomics/*.tsv.gz` - updated Modomics files
- `R/coordinates.R` - major rewrite with auto-classification
- `R/globals.R` - added offset/type globals
- `tests/testthat/test-coordinates.R` - expanded to 46 tests

**Excluded tRNAs:**
- SeC (selenocysteine) - structurally incompatible
- Mitochondrial - different architecture
- Initiator Met - different structural features

**Next steps:**
- Merge feature/global-tRNA-coordinates to devel
- Update vignette with new coordinate groups

---

### 2025-11-24 - Heatmap Type Splitting Enhancement

**Work completed:**
- Split heatmap visualization into Type I and Type II tRNAs
- Type II tRNAs (Leu, Ser, Tyr, SeC) have extended variable loops (9-24 extra nucleotides)
- Plots are stacked vertically using patchwork
- Removed label_interval parameter - now shows ALL Sprinzl labels
- Positions without Sprinzl labels display "NA"
- Added classify_trna_type() helper function to identify Type I vs Type II
- Generated new example plot: man/figures/heatmap-split-example.png
- All 86 tests pass (19 coords + 2 bcerror + 11 plots + 53 structure + 1 nested)

**Files modified:**
- `R/coordinates.R` (added classify_trna_type, lines 180-203)
- `R/plots.R` (refactored plot_bcerror_heatmap, added .plot_heatmap_internal)
- `DESCRIPTION` (added patchwork dependency)
- `tests/testthat/test-plots.R` (new file, 11 tests)
- `man/classify_trna_type.Rd` (new)
- `man/plot_bcerror_heatmap.Rd` (updated)
- `NAMESPACE` (added classify_trna_type export)

**Breaking changes:**
- `label_interval` parameter removed from plot_bcerror_heatmap()
- Return type changed to patchwork when split_by_type=TRUE (default)
- To get old single-plot behavior, use split_by_type=FALSE

**Next steps:**
- Merge to devel
- Update vignette with new plot behavior

---

### 2025-11-24 - Structure Visualization Testing & Finalization

**Work completed:**
- Validated structure visualization with yeast example data
- Generated 5 example plots successfully:
  - Genome-wide heatmap
  - Single tRNA cloverleaf with error rates
  - Aggregated structure (consensus view)
  - Structure with Modomics overlay
  - Faceted comparison (grande vs petite)
- Added 14 comprehensive tests for structure.R functions (53 test assertions total)
- Tests cover: load_structure_template(), load_modomics(), plot_trna_structure()
- Improved roxygen2 documentation for plot_trna_structure():
  - Added cross-reference to plot_bcerror_heatmap()
  - Enhanced parameter descriptions
  - Made examples runnable with clover_example()
  - Added biological interpretation to return value
- Package passes R CMD check: 0 errors, 0 warnings, 2 harmless notes
- All 74 tests pass (19 coordinate + 2 bcerror + 14 structure + 39 nested assertions)

**Files modified:**
- `tests/testthat/test-structure.R` (new - 14 tests covering core functionality)
- `R/structure.R` (enhanced documentation)
- `man/plot_trna_structure.Rd` (regenerated)
- `man/figures/` (4 example PNG plots saved for future documentation)

**Next steps:**
- Merge feature/global-tRNA-coordinates to devel
- Update vignette with complete runnable examples
- Add example plots to README
- Consider Type II tRNA support (already scoped in deferred work)

---

### 2025-11-24 - Structure Visualization Implementation

**Work completed:**
- Implemented `plot_trna_structure()` for tRNA cloverleaf visualization
- Created `R/structure.R` with functions:
  - `plot_trna_structure()` - Cloverleaf secondary structure plot with data overlay
  - `load_structure_template()` - Load consensus Type I tRNA structure coords
  - `load_modomics()` - Load Modomics modification reference data
  - `read_mod_calls()` - Import upstream pipeline consensus modification calls
  - `parse_r2dt_json()` - Internal helper to parse R2DT JSON files
- Bundled Modomics data for E. coli, yeast, and human
- Created consensus structure template from R2DT JSON (Type I tRNA)
- Added support for:
  - Single tRNA or aggregated (consensus) views
  - Modomics overlay showing known modification positions
  - Faceted comparison of multiple conditions
  - Multiple fill options (error_rate, bcerror_residual, p_adjusted, etc.)

**Files added/modified:**
- `R/structure.R` (new)
- `R/globals.R` (added structure viz variables)
- `DESCRIPTION` (added jsonlite, purrr dependencies)
- `inst/extdata/structure/` (consensus template RDS and JSON)
- `inst/extdata/modomics/` (Modomics TSV files for 3 organisms)

**Next steps:**
- Add tests for structure functions
- Add per-tRNA mode using individual R2DT JSONs
- Update vignette with structure visualization examples

---

### 2025-11-24 - Global Coordinate System Integration

**Work completed:**
- Integrated global tRNA coordinate system from tRNAs-in-space project
- Added pre-computed coordinate files for E. coli K12, S. cerevisiae, and H. sapiens
- Created `R/coordinates.R` with functions:
  - `load_global_coords()` - Load pre-computed coordinate mappings
  - `add_global_coords()` - Join bcerror data with global coordinates
  - `get_global_labels()` - Get Sprinzl labels for plotting
  - `get_region_bounds()` - Get structural region boundaries
  - `available_organisms()` - List available organisms
- Added `plot_bcerror_heatmap()` - Heatmap visualization using global coordinates
- Fixed `CloverSE()` class - Now properly stores bcerror and coordinates
- Added comprehensive tests for coordinate functions (19 new tests)
- Package passes R CMD check with 0 errors, 0 warnings

**Files added/modified:**
- `R/coordinates.R` (new)
- `R/plots.R` (added plot_bcerror_heatmap)
- `R/clover-se.R` (rewrote CloverSE constructor)
- `R/globals.R` (added new global variables)
- `inst/extdata/coords/` (new directory with coordinate TSVs)
- `tests/testthat/test-coordinates.R` (new)

**Next steps:**
- Implement `plot_trna_structure()` for secondary structure diagrams
- Implement `calc_diff_exp()` for differential expression
- Implement `calc_diff_mod()` for differential modification
- Update vignette to be runnable
- Consider adding Modomics modification mapping integration

---

### 2024-11-24 - AI Readiness Setup

**Work completed:**
- Added AI assistant configuration files (CLAUDE.md, .cursor/rules.md, .github/copilot-instructions.md)
- Added session workflow instructions
- Created PROJECT_TRACKER.md

**Next steps:**
- Continue package development
- Add more unit tests
- Expand vignette documentation

---

## Known Issues

*No known issues at this time.*

## Deferred Work

### Type II tRNA Support

**Status**: Deferred for future implementation

Type II tRNAs (Leu, Ser, Tyr, SeC) have extended variable arms with 9-24 extra nucleotides that don't fit the standard 76nt consensus template. These tRNAs require:
- Separate Type II consensus template extracted from R2DT JSON
- Extended global coordinate mapping (e-positions)
- Potentially different visualization layout

**Action items for future sessions:**
1. Extract Type II template from an R2DT JSON (e.g., tRNA-Leu)
2. Add `template = "type2"` option to `plot_trna_structure()`
3. Handle e-position mapping in coordinate joins

## Architecture Notes

### Data Model
- Uses `SummarizedExperiment` as base class for `CloverSE`
- Assays store counts and base-calling error data
- Metadata stores file paths (pod5, bam)

### File Organization
- `R/clover-se.R` - Main class definition and constructors
- `R/coordinates.R` - Global coordinate system functions
- `R/plots.R` - Visualization functions
- `R/utils.R` - Helper functions
- `R/read-bcerror.R` - Data input functions
- `R/globals.R` - Global variable declarations
- `inst/extdata/coords/` - Pre-computed global coordinate files

## Useful Commands

```r
# Development workflow
devtools::load_all()
devtools::document()
devtools::test()
devtools::check()

# Build site
pkgdown::build_site()
```
