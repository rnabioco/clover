# Create a SummarizedExperiment from pipeline output.

Read a pipeline configuration file and load all result types into a
[SummarizedExperiment::SummarizedExperiment](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
object. The returned object contains count matrices as assays, tRNA
metadata as rowData, and sample metadata as colData.

## Usage

``` r
create_clover(
  config_path,
  types = c("charging", "bcerror", "odds_ratios"),
  sample_info = NULL,
  min_count = 10
)
```

## Arguments

- config_path:

  Path to a pipeline `config.yaml` file.

- types:

  Character vector of result types to load. Valid values: `"charging"`,
  `"bcerror"`, `"odds_ratios"`.

- sample_info:

  An optional data frame with a `sample_id` column and additional
  experimental factor columns (e.g., `condition`, `replicate`). If
  `NULL`, a minimal colData is created from sample names.

- min_count:

  Minimum total count across all samples for a tRNA to be retained in
  count matrices. Default `10`.

## Value

A
[SummarizedExperiment::SummarizedExperiment](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
with:

- **assay "counts"**: abundance count matrix (charged + uncharged)

- **assay "charging"**: charging count matrix (charged/uncharged columns
  per sample), only if `"charging"` is in `types`

- **colData**: sample metadata

- **metadata**: list with `$config`, `$bcerror`, `$odds_ratios`,
  `$fasta` as available

## Examples

``` r
if (FALSE) { # \dontrun{
se <- create_clover("path/to/config.yaml")
SummarizedExperiment::assay(se, "counts")
SummarizedExperiment::colData(se)
} # }
```
