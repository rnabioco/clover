# Cursor Rules for clover

## Project Context

This is an R package for analyzing nanopore tRNA sequencing data. It uses Bioconductor's SummarizedExperiment framework and tidyverse conventions.

## Code Generation Guidelines

### R Style

- Follow tidyverse style guide
- Use `<-` for assignment, not `=`
- Use snake_case for function and variable names
- Prefer tidyverse verbs (dplyr, tidyr) for data manipulation
- Use ggplot2 for visualization

### Function Documentation

Always include roxygen2 documentation for exported functions:

```r
#' Brief one-line description
#'
#' Longer description if needed.
#'
#' @param param_name Description of parameter
#' @return Description of return value
#' @export
#' @examples
#' # Example usage
#' function_name(arg)
```

### Bioconductor Conventions

- Use S4 classes that extend Bioconductor base classes
- SummarizedExperiment for experiment data
- GRanges/GRangesList for genomic coordinates
- Biostrings for sequence data

### Testing

- Use testthat (edition 3)
- Test files go in `tests/testthat/`
- Name test files `test-*.R`

### Imports

When adding new dependencies:
1. Add to DESCRIPTION Imports or Suggests
2. Use `@import` or `@importFrom` in roxygen for frequently used functions
3. Use `pkg::function()` for occasional use

## File Patterns

- `R/*.R` - Package source code
- `man/*.Rd` - Generated docs (do not edit directly)
- `tests/testthat/test-*.R` - Unit tests
- `vignettes/*.qmd` - Quarto vignettes

## Build Commands

```r
devtools::document()  # Generate documentation
devtools::test()      # Run tests
devtools::check()     # Full R CMD check
```
