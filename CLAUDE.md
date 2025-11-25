# clover

R package for analyzing and visualizing nanopore tRNA sequencing data.

## Project Overview

**Repository**: https://github.com/rnabioco/clover
**Organization**: RNA Bioscience Initiative, University of Colorado Anschutz
**Status**: Active development (v0.0.0.9000)

### Purpose

clover provides tools for:
- Differential expression analysis of tRNAs
- Base-calling error rate analysis and visualization
- tRNA modification detection and comparison
- Visualization including heatmaps and secondary structure diagrams

### Data Model

The package uses `SummarizedExperiment::RangedSummarizedExperiment` to store:
- `counts`: Expression levels per tRNA
- `bcerror`: Base-calling error rates per tRNA and position
- Reference sequences (FASTA)
- Modification annotations

## Development

### Build System

- **Documentation**: roxygen2 with markdown support
- **Testing**: testthat (edition 3)
- **Vignettes**: Quarto and knitr
- **Website**: pkgdown with rbitemplate
- **Formatting**: Posit AIR

### Key Commands

```bash
# Install dependencies
pak::pak("rnabioco/clover")

# Build documentation
devtools::document()

# Run tests
devtools::test()

# Check package
devtools::check()

# Build pkgdown site
pkgdown::build_site()
```

### Branch Strategy

- `devel`: Main development branch (CI target)
- Feature branches for new work

## Code Style

### R Conventions

- Use tidyverse style (dplyr, tidyr, ggplot2)
- roxygen2 for all exported functions
- S4 classes extend Bioconductor base classes
- Pipe-friendly function design

### Documentation

- All exported functions require roxygen2 documentation
- Include `@param`, `@return`, `@export`, and `@examples` tags
- Use markdown in roxygen comments (`Roxygen: list(markdown = TRUE)`)

### Function Patterns

```r
# Exported function with roxygen2
#' Brief description
#'
#' @param arg description
#' @return description
#' @export
#' @examples
#' example_code()
function_name <- function(arg) {
  # implementation
}
```

## Dependencies

### Bioconductor
- SummarizedExperiment
- GenomicRanges
- Biostrings

### Tidyverse
- dplyr, tidyr, readr, ggplot2, forcats

## Testing

Tests live in `tests/testthat/`. Run with:

```r
devtools::test()
testthat::test_file("tests/testthat/test-specific.R")
```

## File Structure

```
R/                  # Source code
man/                # Generated documentation (do not edit)
tests/testthat/     # Unit tests
vignettes/          # Quarto vignettes
inst/extdata/       # Example data files
pkgdown/            # Website configuration
```

## CI/CD

GitHub Actions workflows:
- `R-CMD-check.yaml`: Package checks on Ubuntu (R devel + release)
- `pkgdown.yaml`: Documentation site deployment
