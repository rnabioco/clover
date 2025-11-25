# Project Tracker

Track development progress and session notes for clover.

## Current Status

**Version**: 0.0.0.9000 (active development)
**Branch**: See `git branch` for current working branch

## Recent Sessions

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

## Architecture Notes

### Data Model
- Uses `SummarizedExperiment` as base class for `CloverSE`
- Assays store counts and base-calling error data
- Metadata stores file paths (pod5, bam)

### File Organization
- `R/clover-se.R` - Main class definition and constructors
- `R/plots.R` - Visualization functions
- `R/utils.R` - Helper functions
- `R/read-bcerror.R` - Data input functions

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
