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
config <- read_pipeline_config(clover_example("ecoli/config.yaml"))
list_pipeline_files(config, types = c("charging", "bcerror"))
#> # A tibble: 12 × 3
#>    sample_id    type     path                                                   
#>    <chr>        <chr>    <chr>                                                  
#>  1 wt-15-ctl-01 charging /home/runner/work/_temp/Library/clover/extdata/ecoli/s…
#>  2 wt-15-ctl-02 charging /home/runner/work/_temp/Library/clover/extdata/ecoli/s…
#>  3 wt-15-ctl-03 charging /home/runner/work/_temp/Library/clover/extdata/ecoli/s…
#>  4 wt-15-inf-01 charging /home/runner/work/_temp/Library/clover/extdata/ecoli/s…
#>  5 wt-15-inf-02 charging /home/runner/work/_temp/Library/clover/extdata/ecoli/s…
#>  6 wt-15-inf-03 charging /home/runner/work/_temp/Library/clover/extdata/ecoli/s…
#>  7 wt-15-ctl-01 bcerror  /home/runner/work/_temp/Library/clover/extdata/ecoli/s…
#>  8 wt-15-ctl-02 bcerror  /home/runner/work/_temp/Library/clover/extdata/ecoli/s…
#>  9 wt-15-ctl-03 bcerror  /home/runner/work/_temp/Library/clover/extdata/ecoli/s…
#> 10 wt-15-inf-01 bcerror  /home/runner/work/_temp/Library/clover/extdata/ecoli/s…
#> 11 wt-15-inf-02 bcerror  /home/runner/work/_temp/Library/clover/extdata/ecoli/s…
#> 12 wt-15-inf-03 bcerror  /home/runner/work/_temp/Library/clover/extdata/ecoli/s…
```
