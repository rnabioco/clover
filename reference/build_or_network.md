# Build a network from odds ratio data.

Construct a tidygraph network from pairwise odds ratio or ROR data. Each
node is a position (or isodecoder) and each edge represents a
significant co-occurrence relationship. Requires the tidygraph and
igraph packages.

## Usage

``` r
build_or_network(data, value_col = "ror", min_weight = 0)
```

## Arguments

- data:

  A tibble with columns `pos1`, `pos2`, and a value column specified by
  `value_col`. Optionally includes `p_adj`.

- value_col:

  Column name (string) for edge values. Default `"ror"`.

- min_weight:

  Minimum absolute value of `value_col` for an edge to be included.
  Default `0`.

## Value

A `tbl_graph` object with node metrics `degree` and `betweenness`, or
`NULL` if no edges remain after filtering.

## Examples

``` r
# \donttest{
df <- tibble::tibble(
  pos1 = c("20", "34", "20"),
  pos2 = c("34", "58", "58"),
  ror = c(1.5, -0.8, 0.3)
)
build_or_network(df)
#> # A tbl_graph: 3 nodes and 3 edges
#> #
#> # An undirected simple graph with 1 component
#> #
#> # Node Data: 3 × 3 (active)
#>   position degree betweenness
#>   <chr>     <dbl>       <dbl>
#> 1 20            2           0
#> 2 34            2           0
#> 3 58            2           0
#> #
#> # Edge Data: 3 × 4
#>    from    to edge_type    weight
#>   <int> <int> <chr>         <dbl>
#> 1     1     2 co-occurring    1.5
#> 2     2     3 exclusive       0.8
#> 3     1     3 co-occurring    0.3
# }
```
