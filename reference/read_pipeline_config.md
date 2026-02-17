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
if (FALSE) { # \dontrun{
config <- read_pipeline_config("path/to/config.yaml")
config$samples
config$output_dir
} # }
```
