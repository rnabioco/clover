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

- **colData**: sample metadata

- **metadata**: list with `$config`, `$charging`, `$bcerror`,
  `$odds_ratios`, `$fasta` as available

## Examples

``` r
se <- create_clover(clover_example("ecoli/config.yaml"))
SummarizedExperiment::assay(se, "counts")[1:3, ]
#>                                 wt-15-ctl-01 wt-15-ctl-02 wt-15-ctl-03
#> host-tRNA-Ala-GGC-1-1                    298          397          113
#> host-tRNA-Ala-GGC-1-1-uncharged         4393         5014          883
#> host-tRNA-Ala-GGC-1-2                   1625         1506          265
#>                                 wt-15-inf-01 wt-15-inf-02 wt-15-inf-03
#> host-tRNA-Ala-GGC-1-1                    189          347         1253
#> host-tRNA-Ala-GGC-1-1-uncharged         3640         7133        12837
#> host-tRNA-Ala-GGC-1-2                   1299         1857         4409
SummarizedExperiment::colData(se)
#> DataFrame with 6 rows and 1 column
#>                 sample_id
#>               <character>
#> wt-15-ctl-01 wt-15-ctl-01
#> wt-15-ctl-02 wt-15-ctl-02
#> wt-15-ctl-03 wt-15-ctl-03
#> wt-15-inf-01 wt-15-inf-01
#> wt-15-inf-02 wt-15-inf-02
#> wt-15-inf-03 wt-15-inf-03
```
