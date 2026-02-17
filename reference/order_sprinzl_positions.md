# Order Sprinzl position labels.

Create an ordered factor from Sprinzl position labels. Handles numeric
positions and lettered insertions (e.g., "20a", "20A", "47:e1").

## Usage

``` r
order_sprinzl_positions(labels)
```

## Arguments

- labels:

  Character vector of Sprinzl labels.

## Value

A factor with levels ordered by position number then suffix.

## Examples

``` r
order_sprinzl_positions(c("20a", "20", "1", "21"))
#> [1] 20a 20  1   21 
#> Levels: 1 20 20a 21
```
