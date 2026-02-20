# tRNA structure SVG generation

This directory contains the offline pipeline for generating tRNA cloverleaf
SVG diagrams bundled with the clover package.

## Prerequisites

Install [pixi](https://pixi.sh), then install the `struct` environment:

```bash
pixi install -e struct
```

This provides Python (3.11+), lxml, R2R, and Infernal (cmalign).

## Input files

### FASTA sequences (`fasta/`)

Mature tRNA sequences (without introns) for each organism. These are sourced
from [GtRNAdb](http://gtrnadb.ucsc.edu/) and include only one representative
sequence per gene copy (named `tRNA-{AA}-{Anticodon}-{locus}-{copy}`).

- `eschColi_K_12_MG1655-mature-tRNAs.fa`
- `sacCer3-mature-tRNAs.fa`
- `hg38-mature-tRNAs.fa`

### Covariance models

tRNAscan-SE covariance models for structural alignment with Infernal's
`cmalign`. These models have more consensus columns than the Rfam RF00005
model and properly model the variable arm as a stem in Leu/Ser/Tyr tRNAs.

- `TRNAinf-bact.cm` -- bacterial tRNAs (93 consensus columns)
- `TRNAinf-euk.cm` -- eukaryotic tRNAs (90 consensus columns)
- `TRNAinf-bact-SeC.cm` -- bacterial selenocysteine tRNAs
- `TRNAinf-euk-SeC.cm` -- eukaryotic selenocysteine tRNAs

Source: [tRNAscan-SE](https://github.com/UCSC-LoweLab/tRNAscan-SE) model
library (`lib/models/`).

## Running

From the project root:

```bash
pixi run -e struct python data-raw/structures/generate_svgs.py
```

This will:

1. Align each organism's FASTA against its covariance model using `cmalign`
2. Extract per-sequence secondary structures from the Stockholm alignment
3. Convert to single-sequence Stockholm files with R2R layout directives
4. Run R2R to generate cloverleaf SVG diagrams
5. Parse SVGs to extract nucleotide position metadata (JSON)

Output is written to `inst/extdata/structures/{organism}/`.

## Scripts

- `generate_svgs.py` -- Main pipeline script (alignment, R2R, orchestration)
- `parse_r2r_svg.py` -- SVG parsing and metadata extraction

## Output format

For each tRNA isotype-anticodon pair, two files are generated:

- `{tRNA-name}.svg` -- R2R cloverleaf diagram
- `{tRNA-name}.json` -- Position metadata with fields:
  - `nucleotides`: list of `{pos, base, x, y}` dicts
  - `width`, `height`: SVG dimensions
  - `lines`: base-pair line coordinates
  - `trna_name`: original tRNA identifier
  - `sequence`: RNA sequence
  - `structure`: dot-bracket secondary structure

## Organisms

Organisms are configured in `generate_svgs.py` in the `ORGANISMS` dict. Each
entry specifies the FASTA file and the appropriate covariance model. Only
isotypes with bundled [MODOMICS](https://genesilico.pl/modomics/) data are
included in the output.

| Organism | CM | tRNAs |
|---|---|---|
| *Escherichia coli* K-12 MG1655 | `TRNAinf-bact.cm` | 39 |
| *Saccharomyces cerevisiae* S288c | `TRNAinf-euk.cm` | 38 |
| *Homo sapiens* GRCh38 | `TRNAinf-euk.cm` | 47 |

## Adding a new organism

1. Obtain a mature tRNA FASTA (no introns) and place it in `fasta/`
2. Add an entry to `ORGANISMS` in `generate_svgs.py` with the FASTA path and
   the appropriate CM (`TRNAinf-bact.cm` or `TRNAinf-euk.cm`)
3. Optionally add MODOMICS data to `inst/extdata/modomics/` to filter to
   isotypes with known modifications
4. Re-run the pipeline
