# E. coli test data

Example data from a T4 phage infection time-course experiment on wild-type
*E. coli* K-12 at 15 minutes post-infection. Data were generated with the
[aa-tRNA-seq-pipeline](https://github.com/rnabioco/aa-tRNA-seq-pipeline).

## Samples

6 samples (3 biological replicates per condition):

| Sample | Condition |
|--------|-----------|
| wt-15-ctl-01, wt-15-ctl-02, wt-15-ctl-03 | Uninfected control |
| wt-15-inf-01, wt-15-inf-02, wt-15-inf-03 | T4-infected |

## Files

### Pipeline output (`summary/tables/{sample}/`)

These TSV files are subsets of the full pipeline output, containing 5 host
tRNAs (charged + uncharged = 10 total). They are used by `create_clover()`
to demonstrate the data-loading workflow.

- `{sample}.bcerror.tsv.gz` -- per-position base-calling error rates
- `{sample}.charging.cpm.tsv.gz` -- charging counts per million
- `{sample}.odds_ratios.tsv.gz` -- pairwise modification co-occurrence
- `{sample}.align_stats.tsv.gz` -- alignment statistics

### Preprocessed data

- **`bcerror_summary.rds`** -- Preprocessed bcerror data for 10 tRNAs (5
  most abundant + 5 least abundant), used by the vignette for the
  `plot_mod_heatmap()` and `plot_mod_landscape()` examples. Contains a list
  with two tibbles:
  - `bcerror_charged`: per-sample error rates (ref, pos, error_rate, mis,
    condition, sample_id)
  - `bcerror_summary`: mean error rates across replicates (ref, pos,
    condition, mean_error, mean_mis)

  The 10 tRNAs were chosen to span the abundance range:

  | tRNA | Category |
  |------|----------|
  | host-tRNA-Glu-TTC-1-1 | Most abundant |
  | host-tRNA-Thr-TGT-1-1 | Most abundant |
  | host-tRNA-Ile-GAT-1-1 | Most abundant |
  | host-tRNA-Phe-GAA-1-1 | Most abundant |
  | host-tRNA-Ser-GCT-1-1 | Most abundant |
  | host-tRNA-Gly-CCC-1-1 | Least abundant |
  | host-tRNA-Ala-GGC-1-1 | Least abundant |
  | host-tRNA-Arg-CCG-1-1 | Least abundant |
  | host-tRNA-Arg-CCT-1-1 | Least abundant |
  | host-tRNA-Ser-CGA-1-1 | Least abundant |

### Other files

- `config.yaml` -- pipeline configuration file (used by `create_clover()`)
- `samples.tsv` -- sample sheet
- `validated.fa.gz` -- reference tRNA FASTA sequences

## Provenance

Source data: `/beevol/home/jhessel/devel/rnabioco/2026-phage-infection/results/`

The pipeline TSV subsets were created by filtering the full output to 5
tRNAs. The `bcerror_summary.rds` was created by reading the full bcerror
TSVs for 10 tRNAs, filtering to charged species, and summarizing mean
error rates per position per condition across the 3 replicates.
