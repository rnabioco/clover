# Modification rewiring analysis

``` r
library(clover)
library(dplyr)
```

## Overview

This article demonstrates clover’s modification rewiring analysis
workflow. “Rewiring” refers to changes in modification co-occurrence
patterns between conditions — some position pairs may gain or lose
coordinated modification.

The workflow operates at the **isodecoder** level, collapsing gene
copies to increase statistical power. It uses the same *E. coli* T4
phage infection dataset as the main vignette.

## Load data

``` r
config_path <- clover_example("ecoli/config.yaml")

sample_info <- data.frame(
  sample_id = c(
    "wt-15-ctl-01",
    "wt-15-ctl-02",
    "wt-15-ctl-03",
    "wt-15-inf-01",
    "wt-15-inf-02",
    "wt-15-inf-03"
  ),
  condition = rep(c("ctl", "inf"), each = 3)
)

se <- create_clover(
  config_path,
  types = c("charging", "odds_ratios"),
  sample_info = sample_info
)

or_data <- S4Vectors::metadata(se)$odds_ratios
or_data$condition <- ifelse(grepl("ctl", or_data$sample_id), "ctl", "inf")
```

## Clean and aggregate odds ratios

The raw odds ratios may contain infinite values from zero-cell
contingency tables.
[`clean_odds_ratios()`](https://rnabioco.github.io/clover/reference/clean_odds_ratios.md)
caps these at a reasonable threshold.

``` r
or_clean <- clean_odds_ratios(or_data)
glimpse(or_clean)
#> Rows: 63,823
#> Columns: 10
#> $ ref            <chr> "host-tRNA-Asp-GTC-1-1", "host-tRNA-Asp-GTC-1-1", "host…
#> $ pos1           <dbl> 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20,…
#> $ pos2           <dbl> 23, 26, 28, 31, 32, 35, 36, 37, 39, 40, 43, 45, 47, 48,…
#> $ odds_ratio     <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0…
#> $ log_odds_ratio <dbl> -23.02585, -23.02585, -23.02585, -23.02585, -23.02585, …
#> $ p_value        <dbl> 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1…
#> $ total_obs      <dbl> 1488, 1488, 1488, 1488, 1488, 1488, 1488, 1488, 1488, 1…
#> $ sample_id      <chr> "wt-15-ctl-01", "wt-15-ctl-01", "wt-15-ctl-01", "wt-15-…
#> $ condition      <chr> "ctl", "ctl", "ctl", "ctl", "ctl", "ctl", "ctl", "ctl",…
#> $ log_or_clean   <dbl> -23.02585, -23.02585, -23.02585, -23.02585, -23.02585, …
```

Next, collapse per-gene odds ratios to the isodecoder level. This
averages across gene copies (e.g., `tRNA-Glu-TTC-1-1` and
`tRNA-Glu-TTC-2-1` become `tRNA-Glu-TTC`).

``` r
or_ctl <- or_clean |>
  filter(condition == "ctl") |>
  aggregate_or_isodecoder()

or_inf <- or_clean |>
  filter(condition == "inf") |>
  aggregate_or_isodecoder()

or_ctl
#> # A tibble: 18,030 × 9
#>    isodecoder      pos1  pos2 mean_or mean_log_or sd_log_or min_pval total_reads
#>    <chr>          <dbl> <dbl>   <dbl>       <dbl>     <dbl>    <dbl>       <dbl>
#>  1 host-tRNA-Ala…    18    33       0       -23.0        NA        1         438
#>  2 host-tRNA-Ala…    19    18       0       -23.0        NA        1         438
#>  3 host-tRNA-Ala…    19    33       0       -23.0        NA        1         438
#>  4 host-tRNA-Ala…    19    38       0       -23.0        NA        1         327
#>  5 host-tRNA-Ala…    19    75       0       -23.0        NA        1         327
#>  6 host-tRNA-Ala…    19    76       0       -23.0        NA        1         327
#>  7 host-tRNA-Ala…    19   101       0       -23.0        NA        1         327
#>  8 host-tRNA-Ala…    28    18       0       -23.0        NA        1         438
#>  9 host-tRNA-Ala…    28    19       0       -23.0        NA        1         438
#> 10 host-tRNA-Ala…    28    30       0       -23.0        NA        1         438
#> # ℹ 18,020 more rows
#> # ℹ 1 more variable: n_copies <int>
```

## Compare conditions with relative odds ratios

[`compute_ror_isodecoder()`](https://rnabioco.github.io/clover/reference/compute_ror_isodecoder.md)
computes the relative odds ratio (ROR) between two conditions with
z-score significance testing.

``` r
ror <- compute_ror_isodecoder(or_inf, or_ctl)
ror
#> # A tibble: 13,348 × 15
#>    isodecoder     pos1  pos2 mean_log_or_num se_num mean_log_or_den se_den   ror
#>    <chr>         <dbl> <dbl>           <dbl>  <dbl>           <dbl>  <dbl> <dbl>
#>  1 host-tRNA-Al…    19    33           -23.0     NA           -23.0     NA     0
#>  2 host-tRNA-Al…    28    19           -23.0     NA           -23.0     NA     0
#>  3 host-tRNA-Al…    28    30           -23.0     NA           -23.0     NA     0
#>  4 host-tRNA-Al…    28    31           -23.0     NA           -23.0     NA     0
#>  5 host-tRNA-Al…    28    32           -23.0     NA           -23.0     NA     0
#>  6 host-tRNA-Al…    28    33           -23.0     NA           -23.0     NA     0
#>  7 host-tRNA-Al…    28    36           -23.0     NA           -23.0     NA     0
#>  8 host-tRNA-Al…    28    40           -23.0     NA           -23.0     NA     0
#>  9 host-tRNA-Al…    28    48           -23.0     NA           -23.0     NA     0
#> 10 host-tRNA-Al…    28    54           -23.0     NA           -23.0     NA     0
#> # ℹ 13,338 more rows
#> # ℹ 7 more variables: ror_se <dbl>, z_score <dbl>, p_value <dbl>, p_adj <dbl>,
#> #   ci_lower <dbl>, ci_upper <dbl>, significant <lgl>
```

## Rewiring scores and dimensionality reduction

Build a matrix of ROR values (isodecoders × position pairs) and compute
per-isodecoder rewiring scores.

``` r
mat <- prepare_rewiring_matrix(ror)
scores <- calculate_rewiring_scores(mat)
scores |> arrange(desc(euclidean_magnitude))
#> # A tibble: 4 × 5
#>   isodecoder        euclidean_magnitude mean_abs_change max_abs_change n_nonzero
#>   <chr>                           <dbl>           <dbl>          <dbl>     <dbl>
#> 1 host-tRNA-Arg-ACG                37.4           4.83              10        14
#> 2 host-tRNA-Glu-TTC                22.4           1.78              10         7
#> 3 host-tRNA-Asp-GTC                22.4           1.74              10         6
#> 4 host-tRNA-Gly-GCC                14.1           0.690             10         2
```

Run PCoA for dimensionality reduction and visualize the result.

``` r
pcoa <- perform_pcoa(mat)
plot_pcoa_rewiring(pcoa, scores)
```

![PCoA of tRNA modification rewiring between control and infected
conditions.](rewiring_files/figure-html/fig-pcoa-1.png)

PCoA of tRNA modification rewiring between control and infected
conditions.

## Network visualization

Build a chord diagram from the ROR data to visualize modification
rewiring between positions.

``` r
ror_agg <- ror |>
  dplyr::group_by(pos1, pos2) |>
  dplyr::summarise(log_ror = mean(ror), .groups = "drop")

plot_chord_ror(ror_agg, title = "Modification rewiring: inf vs ctl")
```

![Chord diagram of modification rewiring between control and infected
conditions.](rewiring_files/figure-html/fig-chord-ror-1.png)

Chord diagram of modification rewiring between control and infected
conditions.

## Session info

``` r
sessionInfo()
#> R version 4.5.2 (2025-10-31)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.3 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] dplyr_1.2.0       clover_0.0.0.9000
#> 
#> loaded via a namespace (and not attached):
#>  [1] SummarizedExperiment_1.40.0 shape_1.4.6.1              
#>  [3] circlize_0.4.17             gtable_0.3.6               
#>  [5] xfun_0.56                   bslib_0.10.0               
#>  [7] ggplot2_4.0.2               GlobalOptions_0.1.3        
#>  [9] htmlwidgets_1.6.4           ggrepel_0.9.6              
#> [11] Biobase_2.70.0              lattice_0.22-7             
#> [13] tzdb_0.5.0                  vctrs_0.7.1                
#> [15] tools_4.5.2                 generics_0.1.4             
#> [17] stats4_4.5.2                parallel_4.5.2             
#> [19] tibble_3.3.1                pkgconfig_2.0.3            
#> [21] Matrix_1.7-4                RColorBrewer_1.1-3         
#> [23] S7_0.2.1                    desc_1.4.3                 
#> [25] S4Vectors_0.48.0            lifecycle_1.0.5            
#> [27] stringr_1.6.0               compiler_4.5.2             
#> [29] farver_2.1.2                textshaping_1.0.4          
#> [31] Biostrings_2.78.0           Seqinfo_1.0.0              
#> [33] htmltools_0.5.9             sass_0.4.10                
#> [35] yaml_2.3.12                 pkgdown_2.2.0              
#> [37] pillar_1.11.1               crayon_1.5.3               
#> [39] jquerylib_0.1.4             tidyr_1.3.2                
#> [41] DelayedArray_0.36.0         cachem_1.1.0               
#> [43] abind_1.4-8                 tidyselect_1.2.1           
#> [45] digest_0.6.39               stringi_1.8.7              
#> [47] purrr_1.2.1                 labeling_0.4.3             
#> [49] cowplot_1.2.0               fastmap_1.2.0              
#> [51] grid_4.5.2                  colorspace_2.1-2           
#> [53] cli_3.6.5                   SparseArray_1.10.8         
#> [55] magrittr_2.0.4              S4Arrays_1.10.1            
#> [57] utf8_1.2.6                  readr_2.2.0                
#> [59] withr_3.0.2                 scales_1.4.0               
#> [61] bit64_4.6.0-1               rmarkdown_2.30             
#> [63] XVector_0.50.0              matrixStats_1.5.0          
#> [65] bit_4.6.0                   ragg_1.5.0                 
#> [67] hms_1.1.4                   evaluate_1.0.5             
#> [69] knitr_1.51                  GenomicRanges_1.62.1       
#> [71] IRanges_2.44.0              viridisLite_0.4.3          
#> [73] rlang_1.1.7                 Rcpp_1.1.1                 
#> [75] glue_1.8.0                  BiocGenerics_0.56.0        
#> [77] vroom_1.7.0                 jsonlite_2.0.0             
#> [79] R6_2.6.1                    MatrixGenerics_1.22.0      
#> [81] systemfonts_1.3.1           fs_1.6.6
```
