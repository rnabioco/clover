# Read a bedMethyl file from modkit pileup

Read per-position modification percentages from a bedMethyl file
produced by `modkit pileup`. The 18-column bedMethyl format includes
per-position counts and percentages for each modification type.

## Usage

``` r
read_bedmethyl(path, min_cov = 1)
```

## Arguments

- path:

  Path to a bedMethyl file (may be gzipped).

- min_cov:

  Minimum valid coverage to retain a position. Default `1`.

## Value

A tibble with columns: `ref`, `pos`, `strand`, `mod_code`, `n_valid`,
`percent_mod`, `n_mod`, `n_canonical`.

## Examples

``` r
# read_bedmethyl("sample.bedmethyl.gz")
```
