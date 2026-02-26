# Shorten tRNA name strings for display.

Strips source prefixes (e.g., `host-`, `phage-`), the `tRNA-` prefix,
and the gene-copy suffix (`-\d+-\d+`) to produce compact labels suitable
for plot axes.

## Usage

``` r
shorten_trna_names(x, strip_prefix = "^(host|phage)-")
```

## Arguments

- x:

  Character vector of tRNA name strings.

- strip_prefix:

  Regex pattern for source prefixes to remove. Default
  `"^(host|phage)-"`.

## Value

Character vector of shortened names (e.g., `"host-tRNA-Ser-CGT-1-1"`
becomes `"Ser-CGT"`).

## Examples

``` r
shorten_trna_names("host-tRNA-Ser-CGT-1-1")
#> [1] "Ser-CGT"
shorten_trna_names(c("phage-tRNA-Glu-TTC-2-1", "tRNA-Ala-AGC-3-1"))
#> [1] "Glu-TTC" "Ala-AGC"
```
