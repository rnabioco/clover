# List pipeline output files for each sample.

Given a parsed pipeline configuration, construct expected file paths for
each sample and result type.

## Usage

``` r
list_pipeline_files(
  config,
  types = c("charging", "bcerror", "odds_ratios", "align_stats")
)
```

## Arguments

- config:

  A list returned by
  [`read_pipeline_config()`](https://rnabioco.github.io/clover/reference/read_pipeline_config.md).

- types:

  Character vector of result types. Valid values: `"charging"`,
  `"bcerror"`, `"odds_ratios"`, `"align_stats"`.

## Value

A tibble with columns `sample_id`, `type`, and `path`.

## Examples

``` r
if (FALSE) { # \dontrun{
config <- read_pipeline_config("path/to/config.yaml")
list_pipeline_files(config, types = c("charging", "bcerror"))
} # }
```
