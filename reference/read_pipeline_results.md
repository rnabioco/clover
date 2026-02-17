# Read pipeline results for all samples.

Convenience function that reads a pipeline config, discovers output
files, and loads data for the requested result types. Each result type
is returned as a combined tibble with a `sample_id` column.

## Usage

``` r
read_pipeline_results(
  config_path,
  types = c("charging", "bcerror", "odds_ratios")
)
```

## Arguments

- config_path:

  Path to a Snakemake `config.yaml` file.

- types:

  Character vector of result types to load. Valid values: `"charging"`,
  `"bcerror"`, `"odds_ratios"`.

## Value

A named list of tibbles, one per requested type. Each tibble includes a
`sample_id` column identifying the source sample.

## Examples

``` r
if (FALSE) { # \dontrun{
results <- read_pipeline_results("path/to/config.yaml")
results$charging
results$odds_ratios
} # }
```
