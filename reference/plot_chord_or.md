# Plot a chord diagram of modification co-occurrence.

Display a circlize chord diagram showing significant pairwise
modification co-occurrence (odds ratios) for a single sample. Chords are
colored by direction: positive odds ratios (co-occurring modifications)
vs. negative (mutually exclusive). Chord width reflects the magnitude of
the log odds ratio.

## Usage

``` r
plot_chord_or(
  odds_data,
  or_col = "log_odds_ratio",
  or_cutoff = 0.5,
  p_col = "p_value",
  p_cutoff = 0.05,
  min_obs = 50,
  positive_color = "#D55E00",
  negative_color = "#0072B2",
  sprinzl_coords = NULL,
  title = NULL,
  transparency = 0.4
)
```

## Arguments

- odds_data:

  A tibble of odds ratio data with columns `pos1`, `pos2`,
  `log_odds_ratio` (or column specified by `or_col`), `p_value` (or
  column specified by `p_col`), and `total_obs`.

- or_col:

  Column name (string) for the odds ratio value. Default
  `"log_odds_ratio"`.

- or_cutoff:

  Minimum absolute value of `or_col` for a chord to be drawn. Default
  `0.5`.

- p_col:

  Column name (string) for the p-value. Default `"p_value"`.

- p_cutoff:

  Maximum p-value for a chord to be drawn. Default `0.05`.

- min_obs:

  Minimum number of observations (`total_obs`) for a pair to be
  included. Default `50`.

- positive_color:

  Color for positive odds ratios (co-occurrence). Default `"#D55E00"`
  (vermillion).

- negative_color:

  Color for negative odds ratios (exclusion). Default `"#0072B2"`
  (blue).

- sprinzl_coords:

  An optional tibble from
  [`read_sprinzl_coords()`](https://rnabioco.github.io/clover/reference/read_sprinzl_coords.md)
  used to order sectors by Sprinzl position and color by structural
  region.

- title:

  Optional plot title.

- transparency:

  Transparency for chord colors (0 = opaque, 1 = fully transparent).
  Default `0.4`.

## Value

`invisible(NULL)`. Called for its side effect of producing a base
graphics plot.

## Examples

``` r
if (FALSE) { # \dontrun{
or_data <- read_odds_ratios("sample1.odds_ratios.tsv.gz")
plot_chord_or(or_data)
} # }
```
