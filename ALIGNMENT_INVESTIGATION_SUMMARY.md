# Sprinzl Position Alignment Investigation - Summary of Findings

**Date**: 2025-11-25
**Issue**: Position 55 (and 8/11 common modification sites) misaligned between Type I and Type II tRNAs

## Investigation Summary

### Tests Performed

1. **Current production files**: Confirmed misalignment (Position 55: Type I=87, Type II=98)
2. **Simple `templateNumberingLabel` fix**: Massive collisions (70+ positions)
3. **Dual-system mode**: Creates separate incomparable coordinate spaces
4. **Offset-based grouping**: Still has collisions within offset groups

### Root Cause: Two-Dimensional Problem

The alignment issue stems from **two independent sources of variation** in R2DT output:

#### Dimension 1: Labeling Offset
Different tRNA families have different offsets between `templateNumberingLabel` and `templateResidueIndex`:

**E. coli:**
- Offset -1: 11 tRNAs (Pro, Ser, Tyr families)
- Offset  0: 60 tRNAs (majority)
- Offset +1: 11 tRNAs (various families)

**Yeast:**
- Offset  0: ~50% of tRNAs
- Offset +1: ~50% of tRNAs

#### Dimension 2: Structural Type
- **Type I**: Standard tRNAs (~76 nt)
- **Type II**: Extended variable arm tRNAs (~85-90 nt, includes e1-e24 positions)

### Why Single-Dimension Solutions Fail

**Grouping by offset alone:**
- Offset -1 group contains: Pro (Type I), Ser (Type II), Tyr (Type II)
- These have different structures → collisions between positions 47-76

**Grouping by type alone (dual-system):**
- Within Type I: Multiple offset patterns → different global_index values
- Within Type II: Multiple offset patterns → different global_index values
- Result: Separate coordinate spaces, not aligned

**Neither offset nor label alone:**
- `templateResidueIndex`: Template-specific, misaligns functionally equivalent positions
- `templateNumberingLabel`: Collides due to offset heterogeneity

## Potential Solutions

### Option 1: Multi-Dimensional Grouping (offset × type)

Generate separate coordinate files for each **offset + type** combination:

**E. coli would need:**
- `offset-1_type1.tsv` (Pro)
- `offset-1_type2.tsv` (Ser, Tyr)
- `offset0_type1.tsv` (most Type I tRNAs)
- `offset0_type2.tsv` (most Type II tRNAs)
- `offset+1_type1.tsv`
- `offset+1_type2.tsv`

**Yeast would need:**
- `offset0_type1.tsv`
- `offset0_type2.tsv`
- `offset+1_type1.tsv`
- `offset+1_type2.tsv`

**Pros:**
- Each file would have proper internal alignment
- Using `templateNumberingLabel` within each group should work
- No collisions within groups

**Cons:**
- Requires 4-6 separate coordinate files per organism
- Users must know which file to use for each tRNA
- Cross-group comparisons still not possible
- Complex for downstream analysis (clover would need to handle multiple files)

### Option 2: Offset Normalization

Apply offset correction to standardize all tRNAs to offset=0, then use `templateNumberingLabel`:

```python
# For each tRNA:
# 1. Calculate its offset
# 2. Normalize: adjusted_index = templateResidueIndex + offset
# 3. Use adjusted values as preferred label
```

**Pros:**
- Single coordinate file per organism
- All tRNAs directly comparable
- Simpler for downstream users

**Cons:**
- More complex transformation logic
- Must validate that normalization preserves biological meaning
- Need to verify no new collisions are introduced

### Option 3: Accept Current Limitations

Document that:
- The coordinate system is primarily for within-type comparisons
- Cross-type comparisons at specific positions require manual verification
- Position "55" in Type I ≠ position "55" in Type II in the global coordinate space

**Pros:**
- No code changes needed
- Current system works for many analyses

**Cons:**
- User expectations not met (position 55 doesn't align across types)
- Heatmaps misleading for cross-type comparisons
- Defeats stated purpose of "functionally equivalent positions align"

### Option 4: Investigate R2DT Source

Contact R2DT developers to understand:
- Why do labeling offsets vary?
- Is `templateNumberingLabel` or `templateResidueIndex` the "canonical" coordinate?
- Can R2DT provide a true canonical position that's consistent across templates?

**Pros:**
- Could reveal a proper solution we're missing
- May find R2DT already has this solved

**Cons:**
- External dependency
- May not lead to immediate solution
- R2DT offsets may reflect real structural biology

## Data Points for Decision

### Current Usage Patterns
- clover currently uses single unified coordinate files
- Examples show cross-type heatmaps (but with misalignment)
- Modomics integration expects single coordinate space

### User Impact
- Option 1 (multi-file): Requires clover refactoring to handle multiple coordinate files
- Option 2 (normalization): Transparent to clover users (same file format)
- Option 3 (accept): No changes but limited functionality
- Option 4 (investigate): Unknown timeline

### Technical Feasibility
- Option 1: High (similar to existing dual-system code)
- Option 2: Medium (need to validate normalization logic)
- Option 3: Trivial (documentation only)
- Option 4: Unknown

## Key Discovery: Documentation Error

After careful investigation, **the OUTPUT_FORMAT.md documentation has the field descriptions backwards**:

**Documentation claims:**
- `templateResidueIndex` = canonical/functional
- `templateNumberingLabel` = template-specific alignment

**Actual R2DT behavior:**
- `templateNumberingLabel` = canonical Sprinzl position (e.g., position '18' is always labeled '18')
- `templateResidueIndex` = template-specific index that shifts with insertions (e.g., '18' might be index 18 or 19)

**Evidence:**
- Pro with insertion '17a': label '18' has index=19 (shifted by insertion)
- Ala without insertion: label '18' has index=18 (no shift)
- Both correctly identify canonical position 18 via the LABEL, not the INDEX

## Why Simple Fixes Fail

**Using `templateNumberingLabel` alone**: Still causes collisions despite being canonical
- Reason unclear - possibly due to interpolation artifacts or missing label patterns
- Tested but validation failed with 70+ position collisions

**Offset normalization**: Also causes collisions
- Different tRNAs have incompatible coordinate schemes even after offset correction
- All three offset groups (-1, 0, +1) still collide when using labels

## Recommended Approach: Multi-File System

Since perfect alignment across all tRNAs isn't essential at this phase, implement a practical multi-file system:

### Strategy

Generate **offset × type** coordinate files for each organism:

**E. coli (6 files):**
- `ecoliK12_offset-1_type1.tsv` (~3 Type I tRNAs: Pro family)
- `ecoliK12_offset-1_type2.tsv` (~8 Type II tRNAs: Ser, Tyr families)
- `ecoliK12_offset0_type1.tsv` (~55 Type I tRNAs: majority)
- `ecoliK12_offset0_type2.tsv` (~5 Type II tRNAs: some Leu)
- `ecoliK12_offset+1_type1.tsv` (~8 Type I tRNAs: various families)
- `ecoliK12_offset+1_type2.tsv` (~3 Type II tRNAs: various)

**Yeast (4 files):**
- `sacCer_offset0_type1.tsv`
- `sacCer_offset0_type2.tsv`
- `sacCer_offset+1_type1.tsv`
- `sacCer_offset+1_type2.tsv`

### Implementation in tRNAs-in-space

1. Modify `build_pref_label()` to prefer `templateNumberingLabel` (correct canonical coordinate)
2. Add offset calculation function
3. Add command-line option: `--split-by-offset-and-type`
4. Generate separate files for each combination
5. Each file validated independently (no collisions within groups)

### Integration in clover

**Option A: Metadata-based selection**
- Add offset/type metadata to each tRNA
- Load appropriate coordinate file based on tRNA classification
- Requires clover to detect offset/type for each tRNA

**Option B: Unified file with group column**
- Combine all files with added `coord_group` column
- Users filter by group for analyses
- Simpler integration but larger files

**Option C: Accept current limitations**
- Keep single unified file per organism
- Document cross-type comparison caveats
- Update examples to show within-type analyses

### Trade-offs

**Pros:**
- Each file internally consistent (positions align within group)
- No collisions within offset+type groups
- Uses correct canonical coordinate (`templateNumberingLabel`)
- Pragmatic solution that works today

**Cons:**
- 4-6 files per organism
- Can't directly compare across groups
- More complex for users
- Defeats original "single coordinate space" goal

## Files Generated During Investigation

- `/tmp/test_ecoli_alignment_fix.py` - Test of simple templateNumberingLabel fix
- `/tmp/test_offset_grouping.py` - Test of offset-based grouping
- `/tmp/ecoli_dual_type1.tsv`, `/tmp/ecoli_dual_type2.tsv` - Dual-system test outputs
- This summary document

## User Decision (2025-11-25)

**Confirmed:**
- Perfect alignment NOT essential at this phase
- Multiple coordinate files approach is acceptable
- Wants to see a prototype/plot for E. coli to verify approach works

**Next Steps:**
1. Implement multi-file generation in tRNAs-in-space (offset × type split)
2. Generate E. coli coordinate files
3. Create visualization showing alignment within groups
4. Decide on clover integration approach (Option A, B, or C)
5. Update documentation to clarify `templateNumberingLabel` vs `templateResidueIndex`
