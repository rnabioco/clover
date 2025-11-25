# Sprinzl Position Alignment Issue - Investigation Summary

**Date**: 2024-11-25
**Issue**: Sprinzl position 55 (and other positions) do not align across Type I and Type II tRNAs in heatmaps

## Problem Statement

User noticed that Sprinzl position 55, which has high base-calling error (making it a good "tracer"), does not align vertically in heatmaps when comparing Type I and Type II tRNAs. The high error "signal" appears at different x-axis locations for different tRNA types.

## Investigation Findings

### 1. Current Misalignment Confirmed

Testing the current coordinate files in `inst/extdata/coords/sacCer_global_coords.tsv.gz` shows **8 out of 11 common modification sites are misaligned**:

```
Position 8:  ✓ ALIGNED
Position 13: ✓ ALIGNED
Position 16: ✓ ALIGNED
Position 20: ✗ MISALIGNED (Type I: 21, Type II: 23)
Position 26: ✗ MISALIGNED (Type I: 34, Type II: 36)
Position 34: ✗ MISALIGNED (Type I: 42, Type II: 44)
Position 37: ✗ MISALIGNED (Type I: 45, Type II: 47)
Position 47: ✗ MISALIGNED (Type I: 55, Type II: 72)
Position 54: ✗ MISALIGNED (Type I: 78, Type II: 91)
Position 55: ✗ MISALIGNED (Type I: 79, Type II: 92) ← Reported issue
Position 58: ✗ MISALIGNED (Type I: 82, Type II: 96)
```

Even worse: **Some positions have multiple global_index values within the same tRNA type**. For example, position 55 in Type I tRNAs maps to global_index values 78, 79, 80, AND 81.

### 2. Files Are From tRNAs-in-Space (But Outdated/Wrong)

- Current files in `clover-globalcoords/inst/extdata/coords/` were copied from tRNAs-in-space project
- tRNAs-in-space repository: `/Users/laurawhite/Projects/tRNAs-in-space`
- Documentation claims coordinate system was fixed on Nov 17, 2025 with "zero collisions"
- However, files still show massive misalignment issues

### 3. Root Cause: R2DT Template Index Variations

**Verified from R2DT raw JSON outputs:**

| tRNA | templateNumberingLabel | templateResidueIndex (from R2DT) | global_index (current) |
|------|------------------------|----------------------------------|------------------------|
| nuc-tRNA-Ala-AGC | "55" | **54** | 79 |
| nuc-tRNA-Ala-TGC | "55" | **55** | 80 |
| nuc-tRNA-Leu-CAA | "55" | **66** | 92 |

**Why this happens:**
- R2DT uses different structural templates for different tRNA families
- Templates have insertions/deletions that shift the `templateResidueIndex`
- BUT: `templateNumberingLabel` remains consistent (all are "55")

**Current trnas_in_space.py behavior:**
- `build_pref_label()` function (lines 326-346) prefers `templateResidueIndex` over `templateNumberingLabel`
- Logic: For positions 1-76, use numeric `sprinzl_index` (from `templateResidueIndex`)
- Only uses `sprinzl_label` (from `templateNumberingLabel`) for insertions like "e1", "20a", etc.

```python
# Current logic in trnas_in_space.py
num_ok_str = num_ok.astype("Int64").astype(str).replace({"<NA>": ""})
pref = num_ok_str.mask(num_ok_str.eq(""), other=lbl)  # Prefers numeric index
```

This means:
- Ala-AGC with `templateResidueIndex=54` gets preferred label "54" → ordinal 66 → global_index 79
- Ala-TGC with `templateResidueIndex=55` gets preferred label "55" → ordinal 67 → global_index 80
- Leu-CAA with `templateResidueIndex=66` gets preferred label "66" → ordinal 78 → global_index 92

### 4. Documentation Says This Should Work

From `tRNAs-in-space/docs/OUTPUT_FORMAT.md`:
- Can't use `sprinzl_label` directly due to insertion suffixes (20A, 20B, etc.)
- System uses `sprinzl_ordinal` → `sprinzl_continuous` → `global_index` pipeline
- **Goal**: "Functionally equivalent positions align across all tRNAs"

From `tRNAs-in-space/docs/archive/coordinate-fixes/COORDINATE_SYSTEM_SCOPE.md`:
- Type I and Type II should share **same coordinate range (1-102)**
- Type II extended arms (e1-e24) map to positions 65-85
- Type I has gap at 65-85 for compatibility

From validation report (Nov 17, 2025):
- Claims "zero collisions across ALL organisms"
- S. cerevisiae: 268 tRNAs, **103 positions** ✅
- Status: "PRODUCTION READY"

**BUT**: Current files don't match this description!

### 5. Regeneration Attempts

Regenerated sacCer coordinates using current tRNAs-in-space script:
```bash
cd /Users/laurawhite/Projects/tRNAs-in-space
python scripts/trnas_in_space.py r2dt_outputs/yeast_jsons /tmp/sacCer_global_coords_unified.tsv
```

Output:
- 103 unique positions (matches documentation)
- "No collisions detected" (validation passes)
- **BUT**: Position 55 still has misalignment (79, 80, 92)

**Key finding**: The validation only checks that different `sprinzl_label` values don't collide. It does NOT check that the same `sprinzl_label` gets the same `global_index`!

## Technical Details

### File Locations

**Current (broken) files:**
- `/Users/laurawhite/Projects/clover-globalcoords/inst/extdata/coords/sacCer_global_coords.tsv.gz`
- `/Users/laurawhite/Projects/clover-globalcoords/inst/extdata/coords/ecoliK12_global_coords.tsv.gz`
- `/Users/laurawhite/Projects/clover-globalcoords/inst/extdata/coords/hg38_global_coords.tsv.gz`

**tRNAs-in-space source:**
- Repository: `/Users/laurawhite/Projects/tRNAs-in-space`
- Script: `scripts/trnas_in_space.py`
- R2DT outputs: `r2dt_outputs/{yeast,ecoli,human}_jsons/`
- Current outputs: `outputs/*.tsv` (also have same misalignment)

**Key functions in trnas_in_space.py:**
- `collect_rows_from_json()` (line 157): Extracts R2DT data
- `build_pref_label()` (line 323): Chooses between index vs label - **THIS IS THE PROBLEM**
- `build_global_label_order()` (line 349): Sorts preferred labels
- `make_continuous_for_trna()` (line 355): Creates continuous coordinates
- `validate_no_global_index_collisions()` (line 111): Validates (but misses same-label issue)

## The Core Question

**Should positions with the same `templateNumberingLabel` have the same `global_index`?**

Current behavior: NO - uses template-specific `templateResidueIndex`
Expected for modification analysis: YES - use consistent `templateNumberingLabel`

**Trade-off:**
- Using `templateResidueIndex`: Respects R2DT's template-specific structure
- Using `templateNumberingLabel`: Aligns functionally equivalent positions

From OUTPUT_FORMAT.md:
> "R2DT provides both `sprinzl_index` (canonical/functional) and `sprinzl_label` (alignment-based)"

But the code treats it opposite:
- `sprinzl_index` = `templateResidueIndex` = template-specific (not canonical)
- `sprinzl_label` = `templateNumberingLabel` = canonical functional position

**This seems backwards!**

## Possible Solutions

### Option 1: Fix build_pref_label() Logic
Change priority to prefer `sprinzl_label` for standard positions:
```python
# Instead of: prefer numeric index, fallback to label
# Use: prefer label, only use index for unlabeled positions
```

**Pros**: Fixes alignment issue at source
**Cons**: May break other assumptions in the code

### Option 2: Post-Process Coordinate Files
Create a script that reassigns `global_index` based on `sprinzl_label`:
- All positions with same `sprinzl_label` get same `global_index`
- Preserves original files, creates corrected versions

**Pros**: Non-invasive, easy to test
**Cons**: Band-aid solution, doesn't fix root cause

### Option 3: Verify R2DT Output Understanding
Double-check documentation of what `templateResidueIndex` vs `templateNumberingLabel` mean:
- Are we misunderstanding R2DT's intent?
- Check R2DT documentation for clarification

### Option 4: Contact Original Developer
The tRNAs-in-space validation claims this works correctly, but testing shows it doesn't.
- Possible: Validation was run on different files than what's in outputs/
- Possible: Bug was introduced after validation
- Possible: We're misunderstanding the intended use case

## Current State

**clover-globalcoords:**
- Using broken coordinate files from tRNAs-in-space
- Heatmaps show misaligned positions
- Example plots in `man/figures/` show incorrect alignment
- Tests pass because they don't check cross-type alignment

**tRNAs-in-space:**
- Script exists at `/Users/laurawhite/Projects/tRNAs-in-space/scripts/trnas_in_space.py`
- R2DT outputs exist in `r2dt_outputs/`
- Claims to be "PRODUCTION READY" but has this fundamental issue
- Documentation suggests it should work, but implementation doesn't match

## Next Steps (Options)

1. **Quick investigation**: Check R2DT documentation to understand templateResidueIndex vs templateNumberingLabel
2. **Test fix**: Modify `build_pref_label()` to prefer labels, regenerate, verify alignment
3. **Report issue**: Document this finding and ask for clarification on intended behavior
4. **Workaround**: Create post-processing script for clover that corrects alignment in-memory

## Files Created During Investigation

- `/tmp/sacCer_global_coords_unified.tsv` - Freshly regenerated coordinates (still have issue)
- `scripts/regenerate_example_plots.R` - Script to regenerate plots after fixing coordinates
- Plan file: `/Users/laurawhite/.claude/plans/kind-gliding-graham.md`

## Code Snippet for Verification

Check alignment of common modification sites:
```r
library(dplyr)
devtools::load_all('.')

coords <- load_global_coords('sacCer')
mod_sites <- c("8", "13", "16", "20", "26", "34", "37", "47", "54", "55", "58")

for (pos in mod_sites) {
  pos_data <- coords |>
    filter(sprinzl_label == pos) |>
    mutate(tRNA_type = classify_trna_type(trna_id)) |>
    select(trna_id, tRNA_type, sprinzl_label, global_index) |>
    distinct(tRNA_type, global_index)

  type1_gi <- pos_data$global_index[pos_data$tRNA_type == 'Type I'][1]
  type2_gi <- pos_data$global_index[pos_data$tRNA_type == 'Type II'][1]

  if (!is.na(type1_gi) && !is.na(type2_gi)) {
    if (type1_gi == type2_gi) {
      cat(sprintf("Position %s: ✓ ALIGNED\n", pos))
    } else {
      cat(sprintf("Position %s: ✗ MISALIGNED (Type I: %d, Type II: %d)\n",
                  pos, type1_gi, type2_gi))
    }
  }
}
```
