# tRNA structure SVG generation

This directory contains the offline pipeline for generating tRNA cloverleaf
SVG diagrams bundled with the clover package.

## Prerequisites

- pixi with the `python` environment configured
- R2R (included in the python pixi environment)

## Running

```bash
pixi run -e python python data-raw/structures/generate_svgs.py
```

This will:

1. Download tRNA sequences and secondary structures from gtRNAdb
2. Convert to Stockholm format for R2R
3. Run R2R to generate cloverleaf SVG diagrams
4. Parse SVGs to extract nucleotide position metadata (JSON)

Output is written to `inst/extdata/structures/{organism}/`.

## Files

- `generate_svgs.py` -- Main pipeline script
- `parse_r2r_svg.py` -- SVG parsing and metadata extraction

## Output format

For each tRNA, two files are generated:

- `{tRNA-name}.svg` -- R2R cloverleaf diagram
- `{tRNA-name}.json` -- Position metadata with fields:
  - `nucleotides`: list of `{pos, base, x, y}` dicts
  - `width`, `height`: SVG dimensions
  - `lines`: base-pair line coordinates
  - `trna_name`: original tRNA identifier
  - `sequence`: RNA sequence
  - `structure`: dot-bracket secondary structure

## Organisms

Currently generates SVGs for organisms with bundled MODOMICS data:

- *Escherichia coli* (K-12 MG1655)
- *Saccharomyces cerevisiae* (S288c)
- *Homo sapiens* (GRCh38)
