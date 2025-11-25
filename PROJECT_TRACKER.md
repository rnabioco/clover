# Project Tracker

Track development progress and session notes for clover.

## Current Status

**Version**: 0.0.0.9000 (active development)
**Branch**: See `git branch` for current working branch

## Recent Sessions

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
