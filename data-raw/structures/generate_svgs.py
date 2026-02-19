"""Generate tRNA cloverleaf SVGs using R2R.

Downloads tRNA sequences and secondary structure data from gtRNAdb,
converts to Stockholm format, runs R2R to produce cloverleaf SVGs,
and extracts position metadata as JSON.

Usage:
    pixi run -e python python data-raw/structures/generate_svgs.py

Prerequisites:
    - pixi python environment with lxml
    - R2R available (via pixi python environment)
"""

import json
import re
import subprocess
import sys
import tempfile
import urllib.request
from pathlib import Path

from parse_r2r_svg import parse_r2r_svg, write_metadata

# Project root (two levels up from this script)
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
OUTPUT_DIR = PROJECT_ROOT / "inst" / "extdata" / "structures"
MODOMICS_DIR = PROJECT_ROOT / "inst" / "extdata" / "modomics"

# gtRNAdb organism IDs and their domain classification
ORGANISMS = {
    "Escherichia_coli": {
        "gtrnadb_id": "Esch_coli_K_12_MG1655",
        "domain": "bacteria",
    },
    "Saccharomyces_cerevisiae": {
        "gtrnadb_id": "Scere3",
        "domain": "eukaryota",
    },
    "Homo_sapiens": {
        "gtrnadb_id": "Hsapi38",
        "domain": "eukaryota",
    },
}

# Standard tRNA dot-bracket notation for a canonical cloverleaf
# Used as fallback when we can't get structure from gtRNAdb
CANONICAL_TRNA_STRUCTURE = (
    "(((((((..((((.......)))).(((((.......))))).....(((((.......))))).)))))))...."
)


def main():
    """Generate SVGs for all configured organisms."""
    for org_name, org_info in ORGANISMS.items():
        print(f"\n{'='*60}")
        print(f"Processing {org_name.replace('_', ' ')}")
        print(f"{'='*60}")

        org_dir = OUTPUT_DIR / org_name
        org_dir.mkdir(parents=True, exist_ok=True)

        # Step 1: Get tRNA sequences and structures from gtRNAdb
        trnas = download_gtrnadb(org_info["gtrnadb_id"], org_info["domain"])

        if not trnas:
            print(f"  WARNING: No tRNA data retrieved for {org_name}")
            continue

        # Step 2: Filter to isotypes with MODOMICS data
        modomics_isotypes = get_modomics_isotypes(org_name)
        if modomics_isotypes:
            filtered = {
                k: v for k, v in trnas.items()
                if extract_isotype(k) in modomics_isotypes
            }
            print(
                f"  Filtered to {len(filtered)}/{len(trnas)} tRNAs "
                f"matching MODOMICS isotypes"
            )
            trnas = filtered if filtered else trnas

        # Step 3: Generate SVGs for each tRNA
        for trna_name, trna_data in trnas.items():
            safe_name = sanitize_filename(trna_name)
            svg_path = org_dir / f"{safe_name}.svg"
            json_path = org_dir / f"{safe_name}.json"

            if svg_path.exists() and json_path.exists():
                print(f"  Skipping {safe_name} (already exists)")
                continue

            print(f"  Generating {safe_name}...")

            try:
                # Create Stockholm file and run R2R
                svg_content = run_r2r(
                    trna_name,
                    trna_data["sequence"],
                    trna_data["structure"],
                )

                if svg_content:
                    # Write SVG
                    svg_path.write_text(svg_content)

                    # Parse SVG and write metadata
                    metadata = parse_r2r_svg(svg_path)
                    metadata["trna_name"] = trna_name
                    metadata["sequence"] = trna_data["sequence"]
                    metadata["structure"] = trna_data["structure"]
                    write_metadata(metadata, json_path)
                    print(f"    OK ({len(metadata['nucleotides'])} nucleotides)")
                else:
                    print(f"    FAILED: R2R produced no output")

            except Exception as e:
                print(f"    ERROR: {e}")

    print(f"\nDone. Output in {OUTPUT_DIR}")


def download_gtrnadb(gtrnadb_id: str, domain: str) -> dict:
    """Download tRNA data from gtRNAdb.

    Tries the current gtRNAdb URLs. Falls back to parsing the
    tRNAscan-SE secondary structure file (.ss format).

    Returns dict mapping tRNA names to {"sequence": ..., "structure": ...}.
    """
    base_url = f"https://gtrnadb.ucsc.edu/genomes/{domain}/{gtrnadb_id}"

    # Try to get the .ss file (secondary structure)
    ss_url = f"{base_url}/{gtrnadb_id}-tRNAs.ss"
    print(f"  Fetching {ss_url}")

    try:
        req = urllib.request.Request(
            ss_url,
            headers={"User-Agent": "clover/0.1 (R package)"},
        )
        with urllib.request.urlopen(req, timeout=30) as resp:
            ss_data = resp.read().decode("utf-8")
        return parse_trnascan_ss(ss_data)
    except Exception as e:
        print(f"  Could not fetch .ss file: {e}")

    # Try the FASTA file as fallback
    fa_url = f"{base_url}/{gtrnadb_id}-mature-tRNAs.fa"
    print(f"  Trying FASTA fallback: {fa_url}")

    try:
        req = urllib.request.Request(
            fa_url,
            headers={"User-Agent": "clover/0.1 (R package)"},
        )
        with urllib.request.urlopen(req, timeout=30) as resp:
            fa_data = resp.read().decode("utf-8")
        return parse_fasta_with_default_structure(fa_data)
    except Exception as e:
        print(f"  Could not fetch FASTA file: {e}")

    return {}


def parse_trnascan_ss(ss_data: str) -> dict:
    """Parse tRNAscan-SE .ss output format.

    Each entry looks like:
        name.trna1 (start-end)
        Length: 73 bp
        Type: Ala
        Anticodon: AGC at 34-36 (34-36)
        Score: 74.48
                 *    |    *    |    *    |
        Seq: GGGCGTGTGGCGTAGTCGG...
        Str: >>>>>.>..>>>>........<<<<.>>>>>...

    The Str line uses > for 5' base pairs and < for 3' base pairs.
    """
    trnas = {}
    entries = re.split(r"\n(?=\S+\.trna\d+\s)", ss_data)

    for entry in entries:
        if not entry.strip():
            continue

        # Extract name
        name_match = re.match(r"(\S+)", entry)
        if not name_match:
            continue

        # Extract type and anticodon
        type_match = re.search(r"Type:\s+(\S+)", entry)
        ac_match = re.search(r"Anticodon:\s+(\S+)", entry)

        if not type_match or not ac_match:
            continue

        aa_type = type_match.group(1)
        anticodon = ac_match.group(1)
        trna_name = f"tRNA-{aa_type}-{anticodon}"

        # Extract sequence
        seq_match = re.search(r"Seq:\s+(\S+)", entry)
        str_match = re.search(r"Str:\s+(\S+)", entry)

        if not seq_match or not str_match:
            continue

        sequence = seq_match.group(1).upper().replace("T", "U")
        structure_raw = str_match.group(1)

        # Convert tRNAscan-SE bracket notation to dot-bracket
        structure = convert_trnascan_to_dotbracket(structure_raw)

        # Keep first occurrence of each isotype-anticodon combo
        if trna_name not in trnas:
            trnas[trna_name] = {
                "sequence": sequence,
                "structure": structure,
            }

    return trnas


def convert_trnascan_to_dotbracket(raw: str) -> str:
    """Convert tRNAscan-SE structure notation to dot-bracket.

    tRNAscan uses > for 5' strand of stems and < for 3' strand.
    We convert to standard ( and ) notation.
    """
    result = []
    for ch in raw:
        if ch == ">":
            result.append("(")
        elif ch == "<":
            result.append(")")
        elif ch == ".":
            result.append(".")
        else:
            result.append(".")
    return "".join(result)


def parse_fasta_with_default_structure(fa_data: str) -> dict:
    """Parse FASTA and assign canonical tRNA structure."""
    trnas = {}
    current_name = None
    current_seq = []

    for line in fa_data.strip().split("\n"):
        if line.startswith(">"):
            if current_name and current_seq:
                seq = "".join(current_seq).upper().replace("T", "U")
                name = extract_trna_name(current_name)
                if name and name not in trnas:
                    structure = fit_canonical_structure(seq)
                    trnas[name] = {
                        "sequence": seq,
                        "structure": structure,
                    }
            current_name = line[1:].strip()
            current_seq = []
        else:
            current_seq.append(line.strip())

    # Last entry
    if current_name and current_seq:
        seq = "".join(current_seq).upper().replace("T", "U")
        name = extract_trna_name(current_name)
        if name and name not in trnas:
            structure = fit_canonical_structure(seq)
            trnas[name] = {
                "sequence": seq,
                "structure": structure,
            }

    return trnas


def extract_trna_name(header: str) -> str | None:
    """Extract a tRNA-AA-Anticodon name from a FASTA header."""
    match = re.search(r"(tRNA-\w+-\w+)", header)
    return match.group(1) if match else None


def extract_isotype(trna_name: str) -> str | None:
    """Extract amino acid isotype from tRNA name."""
    match = re.match(r"tRNA-(\w+)-", trna_name)
    return match.group(1) if match else None


def fit_canonical_structure(sequence: str) -> str:
    """Fit canonical tRNA secondary structure to a sequence.

    Adjusts the canonical dot-bracket structure to match the
    sequence length, primarily by modifying the variable loop region.
    """
    n = len(sequence)
    canonical = CANONICAL_TRNA_STRUCTURE

    if n == len(canonical):
        return canonical
    elif n < len(canonical):
        # Truncate from the variable loop region (positions ~44-48)
        diff = len(canonical) - n
        # Remove from the variable loop (between anticodon arm and T-arm)
        vl_start = 44
        return canonical[:vl_start] + canonical[vl_start + diff:]
    else:
        # Extend the variable loop
        diff = n - len(canonical)
        vl_start = 44
        return canonical[:vl_start] + "." * diff + canonical[vl_start:]


def get_modomics_isotypes(org_name: str) -> set | None:
    """Get set of tRNA isotypes with MODOMICS modification data.

    Reads the cached MODOMICS RDS files to determine which amino acid
    types have known modifications. Returns None if no cached data.
    """
    rds_path = MODOMICS_DIR / f"{org_name}.rds"
    if not rds_path.exists():
        return None

    # We can't read RDS from Python directly, so we extract isotype
    # info by running a quick R command
    try:
        result = subprocess.run(
            [
                "Rscript", "-e",
                f'x <- readRDS("{rds_path}"); '
                f'cat(unique(x$subtype), sep="\\n")',
            ],
            capture_output=True,
            text=True,
            timeout=30,
        )
        if result.returncode == 0 and result.stdout.strip():
            return set(result.stdout.strip().split("\n"))
    except Exception:
        pass

    return None


def sanitize_filename(name: str) -> str:
    """Create a filesystem-safe filename from a tRNA name."""
    return re.sub(r"[^a-zA-Z0-9_-]", "_", name)


def run_r2r(trna_name: str, sequence: str, structure: str) -> str | None:
    """Run R2R on a single tRNA to produce an SVG.

    Creates Stockholm and meta files in a temp directory, runs R2R,
    and returns the SVG content as a string.
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)

        # Write Stockholm file
        sto_path = tmpdir / "trna.sto"
        sto_content = create_stockholm(trna_name, sequence, structure)
        sto_path.write_text(sto_content)

        # Write R2R meta file
        meta_path = tmpdir / "trna.r2r_meta"
        meta_path.write_text(f"{sto_path}\n")

        # Output SVG path
        svg_path = tmpdir / "trna.svg"

        # Run R2R
        try:
            result = subprocess.run(
                ["r2r", "--disable-usage-warning", str(meta_path), str(svg_path)],
                capture_output=True,
                text=True,
                timeout=60,
            )

            if result.returncode != 0:
                print(f"    R2R stderr: {result.stderr[:200]}")
                return None

            if svg_path.exists():
                return svg_path.read_text()
            return None

        except FileNotFoundError:
            print("    ERROR: r2r not found. Install via: pixi install -e python")
            return None
        except subprocess.TimeoutExpired:
            print("    ERROR: R2R timed out")
            return None


def create_stockholm(name: str, sequence: str, structure: str) -> str:
    """Create a Stockholm alignment file for R2R.

    R2R expects standard Stockholm format with SS_cons annotation.
    For single molecules, we provide one sequence with matching structure.
    """
    # Ensure sequence and structure are same length
    seq_len = len(sequence)
    str_len = len(structure)

    if str_len < seq_len:
        structure = structure + "." * (seq_len - str_len)
    elif str_len > seq_len:
        structure = structure[:seq_len]

    # Replace U with T for DNA-style Stockholm (R2R convention)
    seq_display = sequence.replace("U", "T").replace("u", "t")

    # Use a safe sequence ID (no spaces)
    seq_id = sanitize_filename(name)

    lines = [
        "# STOCKHOLM 1.0",
        "",
        f"{seq_id}      {seq_display}",
        f"#=GC SS_cons    {structure}",
        "",
        "//",
    ]
    return "\n".join(lines)


if __name__ == "__main__":
    main()
