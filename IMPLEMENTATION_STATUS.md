# Heatmap Split Implementation - Status Document

**Date**: 2025-11-24
**Branch**: feature/global-tRNA-coordinates
**Task**: Split heatmap by tRNA Type I/II with complete Sprinzl labeling

## What Has Been Completed ✅

### 1. Helper Function - `classify_trna_type()`
**File**: `R/coordinates.R` (lines 180-203)
- Added and exported new function
- Identifies Type II tRNAs: Leu, Ser, Tyr, SeC (extended variable loop)
- All others are Type I (standard)
- Pattern matching: `grepl("tRNA-(Leu|Ser|Tyr|SeC)", trna_ids)`

### 2. Internal Plotting Function - `.plot_heatmap_internal()`
**File**: `R/plots.R` (lines 3-93)
- Extracted core plotting logic to avoid duplication
- **KEY CHANGE**: Shows ALL Sprinzl labels (removed label_interval)
- **KEY CHANGE**: Converts "-1" → "NA" for missing Sprinzl labels
- Accepts optional title_suffix parameter
- Marked @noRd (internal use only)

### 3. Main Function Refactor - `plot_bcerror_heatmap()`
**File**: `R/plots.R` (lines 96-166)

**Signature Change**:
```r
# OLD
plot_bcerror_heatmap(bcerror, value, show_regions, label_interval = 5)

# NEW
plot_bcerror_heatmap(bcerror, value, show_regions, split_by_type = TRUE)
```

**Behavior**:
- `split_by_type = TRUE` (default): Returns patchwork object with Type I and Type II stacked vertically
- `split_by_type = FALSE`: Returns single ggplot object with all tRNAs

**Implementation**:
- Adds `trna_type` column via `classify_trna_type()`
- Filters into type1_data and type2_data
- Calls `.plot_heatmap_internal()` for each
- Uses patchwork `p_type1 / p_type2` to stack

### 4. Dependencies
**File**: `DESCRIPTION`
- Added `patchwork` to Imports (line 24)

### 5. Documentation
**Files**: `man/classify_trna_type.Rd`, `man/plot_bcerror_heatmap.Rd`
- Generated via devtools::document()
- Updated roxygen2 for both functions
- Examples show new usage patterns

### 6. Tests
**File**: `tests/testthat/test-plots.R` (NEW FILE)
- 11 tests total, all passing ✅
- Tests `classify_trna_type()` identification
- Tests split vs non-split behavior
- Tests patchwork return type
- Tests axis labeling (>50 positions shown)
- Tests error handling for missing coordinates

**Test Results**: 86 tests pass total (19 coords + 2 bcerror + 11 plots + 53 structure + 1 nested)

## What Still Needs To Be Done ⚠️

### 7. Generate Example Plots
**Task**: Create new example plots showing the split behavior
**Files to create**:
- `man/figures/heatmap-split-example.png` - Stacked Type I/II plot
- Optionally replace old `heatmap-example.png`

**Commands**:
```r
library(devtools)
load_all(".")
library(dplyr)

bcerr <- read_bcerror(clover_example("yeast/grande.bcerr.tsv.gz")) |>
  add_global_coords("sacCer") |>
  filter(grepl("^nuc-", ref))

# Generate stacked plot
p <- plot_bcerror_heatmap(bcerr, value = "error_rate", show_regions = TRUE)

# Save
ggsave("man/figures/heatmap-split-example.png", p, width = 12, height = 10, dpi = 150, bg = "white")
```

### 8. Run R CMD Check
**Command**: `Rscript -e "devtools::check(vignettes = FALSE)"`
**Expected**: 0 errors, 0 warnings (notes OK)

### 9. Update PROJECT_TRACKER.md
**Section to add**:
```markdown
### 2025-11-24 - Heatmap Type Splitting Enhancement

**Work completed:**
- Split heatmap visualization into Type I and Type II tRNAs
- Type II tRNAs (Leu, Ser, Tyr, SeC) have extended variable loops
- Plots are stacked vertically using patchwork
- Removed label_interval parameter - now shows ALL Sprinzl labels
- Positions without Sprinzl labels display "NA"
- Added classify_trna_type() helper function
- All 86 tests pass

**Files modified:**
- R/coordinates.R (added classify_trna_type)
- R/plots.R (refactored plot_bcerror_heatmap)
- DESCRIPTION (added patchwork dependency)
- tests/testthat/test-plots.R (new, 11 tests)

**Breaking changes:**
- label_interval parameter removed
- Return type changed to patchwork when split_by_type=TRUE

**Next steps:**
- Merge to devel
- Update vignette with new plot behavior
```

### 10. Commit and Push
**Commands**:
```bash
git add R/coordinates.R R/plots.R DESCRIPTION tests/testthat/test-plots.R man/ NAMESPACE PROJECT_TRACKER.md man/figures/
git commit -m "Split heatmap by tRNA type with complete labeling

- Add classify_trna_type() to identify Type I vs Type II tRNAs
- Refactor plot_bcerror_heatmap() to return stacked plots
- Type II tRNAs (Leu, Ser, Tyr, SeC) shown separately
- Remove label_interval parameter - show ALL Sprinzl labels
- Convert -1 to NA for missing labels
- Add patchwork dependency for stacking
- Add 11 new tests, all passing (86 total)

BREAKING CHANGES:
- label_interval parameter removed
- Returns patchwork object by default (set split_by_type=FALSE for old behavior)

🤖 Generated with [Claude Code](https://claude.com/claude-code)

Co-Authored-By: Claude <noreply@anthropic.com>"

git push origin feature/global-tRNA-coordinates
```

## Key Technical Details

### Type I vs Type II Identification
**Type II tRNAs** have extended variable loops (9-24 extra nucleotides):
- Leucine (Leu)
- Serine (Ser)
- Tyrosine (Tyr)
- Selenocysteine (SeC)

**Type I tRNAs**: All others (Ala, Arg, Asn, Asp, Cys, Gln, Glu, Gly, His, Ile, Lys, Met, Phe, Pro, Thr, Trp, Val)

### Axis Labeling Logic
**Old behavior** (removed):
- Showed every 5th Sprinzl label
- Filtered out "-1" positions entirely

**New behavior**:
- Shows EVERY global_index position
- Converts "-1" to "NA" string
- Creates full label vector for all indices

**Implementation** (R/plots.R lines 17-31):
```r
all_indices <- sort(unique(plot_data$global_index))

labels_lookup <- plot_data |>
  dplyr::distinct(global_index, sprinzl_label) |>
  dplyr::mutate(
    sprinzl_label = ifelse(sprinzl_label == "-1", "NA", sprinzl_label)
  )

axis_labels <- sapply(all_indices, function(idx) {
  label <- labels_lookup$sprinzl_label[labels_lookup$global_index == idx]
  if (length(label) == 0) "NA" else label[1]
})
names(axis_labels) <- all_indices
```

### Stacking with Patchwork
**Implementation** (R/plots.R line 161):
```r
p_type1 / p_type2  # Vertical stacking operator
```

Returns a patchwork object that can be printed directly or further customized:
```r
library(patchwork)
p <- plot_bcerror_heatmap(bcerr)
p  # Shows stacked plot
p & theme(axis.text.x = element_text(size = 5))  # Apply theme to both
```

## Backward Compatibility Notes

### Breaking Changes
1. **Parameter removed**: `label_interval` - any code using this will error
2. **Return type changed**: Now returns patchwork (when split_by_type=TRUE) instead of ggplot

### Migration Guide
**Old code**:
```r
p <- plot_bcerror_heatmap(bcerr, label_interval = 10)
print(p)
```

**New code**:
```r
# For split behavior (new default)
p <- plot_bcerror_heatmap(bcerr)  # split_by_type=TRUE implicit
print(p)  # Shows stacked plots

# For old single-plot behavior
p <- plot_bcerror_heatmap(bcerr, split_by_type = FALSE)
print(p)
```

## Files Modified

### Source Files
- **R/coordinates.R** - Added classify_trna_type() (lines 180-203)
- **R/plots.R** - Refactored plot_bcerror_heatmap(), added .plot_heatmap_internal() (lines 3-166)
- **DESCRIPTION** - Added patchwork to Imports (line 24)

### Documentation
- **man/classify_trna_type.Rd** - New (auto-generated)
- **man/plot_bcerror_heatmap.Rd** - Updated (auto-generated)
- **NAMESPACE** - Updated with classify_trna_type export (auto-generated)

### Tests
- **tests/testthat/test-plots.R** - New file, 11 tests

### To Update
- **PROJECT_TRACKER.md** - Add session notes
- **man/figures/** - Add new example plots

## How to Resume Work

1. **Check current status**:
   ```bash
   git status
   devtools::test()  # Should show 86 passing tests
   ```

2. **Generate example plots** (see section 7 above)

3. **Run R CMD check**:
   ```bash
   Rscript -e "devtools::check(vignettes = FALSE)"
   ```

4. **Update PROJECT_TRACKER.md** (see section 9 above)

5. **Commit and push** (see section 10 above)

6. **Merge to devel** (from main worktree):
   ```bash
   cd /Users/laurawhite/Projects/clover  # Main worktree
   git checkout devel
   git merge --no-ff feature/global-tRNA-coordinates
   git push origin devel
   ```

## Current Git State

**Branch**: feature/global-tRNA-coordinates
**Uncommitted changes**:
- R/coordinates.R (classify_trna_type added)
- R/plots.R (refactored)
- DESCRIPTION (patchwork added)
- tests/testthat/test-plots.R (new)
- man/*.Rd (auto-generated, not committed yet)
- NAMESPACE (auto-generated, not committed yet)

## Testing Quick Reference

**Run all tests**:
```r
devtools::test()
```

**Run plot tests only**:
```r
testthat::test_file("tests/testthat/test-plots.R")
```

**Interactively test the function**:
```r
devtools::load_all()
bcerr <- read_bcerror(clover_example("yeast/grande.bcerr.tsv.gz")) |>
  add_global_coords("sacCer")

# Stacked (new default)
p1 <- plot_bcerror_heatmap(bcerr)
print(p1)

# Single plot (old behavior)
p2 <- plot_bcerror_heatmap(bcerr, split_by_type = FALSE)
print(p2)
```
