# Plot a chord diagram of modification rewiring between conditions.

Display a circlize chord diagram showing position pairs where
modification co-occurrence has changed between two conditions, based on
the ratio of odds ratios (ROR). Chords are colored by direction: gained
dependencies (positive log ROR) vs. lost dependencies (negative log
ROR).

## Usage

``` r
plot_chord_ror(
  ror_data,
  ror_cutoff = 0.5,
  gained_color = "#D55E00",
  lost_color = "#0072B2",
  sprinzl_coords = NULL,
  mods = NULL,
  title = NULL,
  transparency = 0.4
)
```

## Arguments

- ror_data:

  A tibble from
  [`compute_ror()`](https://rnabioco.github.io/clover/reference/compute_ror.md)
  with columns `pos1`, `pos2`, and `log_ror`.

- ror_cutoff:

  Minimum absolute value of `log_ror` for a chord to be drawn. Default
  `0.5`.

- gained_color:

  Color for gained dependencies (positive log ROR). Default `"#D55E00"`
  (vermillion).

- lost_color:

  Color for lost dependencies (negative log ROR). Default `"#0072B2"`
  (blue).

- sprinzl_coords:

  An optional tibble from
  [`read_sprinzl_coords()`](https://rnabioco.github.io/clover/reference/read_sprinzl_coords.md)
  used to order sectors and color by structural region.

- mods:

  An optional tibble of modification annotations (e.g., from
  [`fetch_modomics_mods()`](https://rnabioco.github.io/clover/reference/fetch_modomics_mods.md))
  with columns `pos` and `mod1`. When provided along with
  `sprinzl_coords`, modification positions are highlighted as an
  annotation ring.

- title:

  Optional plot title.

- transparency:

  Transparency for chord colors. Default `0.4`.

## Value

`invisible(NULL)`. Called for its side effect of producing a base
graphics plot.

## Examples

``` r
if (FALSE) { # \dontrun{
ror <- compute_ror(combined_or, numerator = "mut", denominator = "wt")
plot_chord_ror(ror)
} # }
```
