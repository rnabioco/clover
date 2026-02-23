# Read a Snakemake pipeline configuration file.

Parse a YAML config file from the tRNA sequencing pipeline and return a
structured list with sample information, output directory, and FASTA
reference path. Relative paths are resolved relative to the config file
directory.

## Usage

``` r
read_pipeline_config(config_path)
```

## Arguments

- config_path:

  Path to a Snakemake `config.yaml` file.

## Value

A list with elements:

- `samples`: A tibble with columns `sample_id` and `data_path`.

- `output_dir`: Resolved path to pipeline output directory.

- `fasta`: Resolved path to reference FASTA file.

- `config_dir`: Directory containing the config file.

## Examples

``` r
config <- read_pipeline_config(clover_example("ecoli/config.yaml"))
config$samples
#> # A tibble: 6 × 2
#>   sample_id    data_path         
#>   <chr>        <chr>             
#> 1 wt-15-ctl-01 /data/wt-15-ctl-01
#> 2 wt-15-ctl-02 /data/wt-15-ctl-02
#> 3 wt-15-ctl-03 /data/wt-15-ctl-03
#> 4 wt-15-inf-01 /data/wt-15-inf-01
#> 5 wt-15-inf-02 /data/wt-15-inf-02
#> 6 wt-15-inf-03 /data/wt-15-inf-03
config$output_dir
#> [1] "/home/runner/work/_temp/Library/clover/extdata/ecoli"
```
