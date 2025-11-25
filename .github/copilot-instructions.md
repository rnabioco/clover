# GitHub Copilot Instructions for clover

## Project Overview

clover is an R/Bioconductor package for analyzing nanopore tRNA sequencing data.

## Code Style

- Use tidyverse style guide conventions
- Assignment with `<-`, not `=`
- snake_case for functions and variables
- Document all exported functions with roxygen2

## Key Patterns

### roxygen2 Documentation

```r
#' Brief description
#'
#' @param arg Description
#' @return Description
#' @export
#' @examples
#' example()
```

### Data Manipulation

Use tidyverse verbs:
```r
data |>
  dplyr::filter(condition) |>
  dplyr::mutate(new_col = transformation) |>
  dplyr::select(col1, col2)
```

### Visualization

Use ggplot2:
```r
ggplot(data, aes(x = x, y = y)) +
  geom_point() +
  theme_minimal()
```

### Bioconductor Classes

- SummarizedExperiment for experiment data
- GRanges for genomic coordinates
- Biostrings for DNA/RNA sequences

## Testing

Use testthat edition 3:
```r
test_that("description", {
  expect_equal(actual, expected)
})
```

## Dependencies

Core: SummarizedExperiment, GenomicRanges, Biostrings, dplyr, tidyr, ggplot2
