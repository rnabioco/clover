"""Generate tRNA cloverleaf SVGs using R2R.

Extracts tRNA sequences and secondary structure data from local gtRNAdb
archive files, converts to Stockholm format, runs R2R to produce
cloverleaf SVGs, and extracts position metadata as JSON.

Usage:
    python data-raw/structures/generate_svgs.py

Prerequisites:
    - Python with lxml
    - R2R on PATH
    - Local gtRNAdb archive tarballs in ~/
"""

import re
import subprocess
import tarfile
import tempfile
from pathlib import Path

from parse_r2r_svg import parse_r2r_svg, write_metadata

# Project root (two levels up from this script)
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
OUTPUT_DIR = PROJECT_ROOT / "inst" / "extdata" / "structures"
MODOMICS_DIR = PROJECT_ROOT / "inst" / "extdata" / "modomics"

HOME = Path.home()

# r2r binary — multistem_junction_bulgey doesn't need the NLOPT solver,
# so we can use whichever r2r is on PATH
R2R_BIN = "r2r"

# Local archive files and the .ss file within each
ORGANISMS = {
    "Escherichia_coli": {
        "archive": HOME / "eschColi_K_12_MG1655-tRNAs.tar.gz",
        "ss_member": "eschColi_K_12_MG1655-tRNAs.ss.sort",
    },
    "Saccharomyces_cerevisiae": {
        "archive": HOME / "sacCer3-tRNAs.tar.gz",
        "ss_member": "sacCer3-tRNAs.ss.sort",
    },
    "Homo_sapiens": {
        "archive": HOME / "hg38-tRNAs.tar.gz",
        "ss_member": "hg38-tRNAs-confidence-set.ss",
    },
}


def main():
    """Generate SVGs for all configured organisms."""
    for org_name, org_info in ORGANISMS.items():
        print(f"\n{'='*60}")
        print(f"Processing {org_name.replace('_', ' ')}")
        print(f"{'='*60}")

        org_dir = OUTPUT_DIR / org_name
        org_dir.mkdir(parents=True, exist_ok=True)

        # Step 1: Extract tRNA sequences and structures from local archive
        trnas = extract_trnas_from_archive(
            org_info["archive"], org_info["ss_member"]
        )

        if not trnas:
            print(f"  WARNING: No tRNA data extracted for {org_name}")
            continue

        print(f"  Found {len(trnas)} unique tRNA isotype-anticodon combos")

        # Step 2: Clean output directory of stale files
        for old_file in org_dir.glob("*.svg"):
            old_file.unlink()
        for old_file in org_dir.glob("*.json"):
            old_file.unlink()

        # Step 3: Filter to isotypes with MODOMICS data
        modomics_isotypes = get_modomics_isotypes(org_name)
        if modomics_isotypes:
            filtered = {
                k: v
                for k, v in trnas.items()
                if extract_isotype(k) in modomics_isotypes
            }
            print(
                f"  Filtered to {len(filtered)}/{len(trnas)} tRNAs "
                f"matching MODOMICS isotypes"
            )
            trnas = filtered if filtered else trnas

        # Step 4: Generate SVGs for each tRNA
        for trna_name, trna_copies in trnas.items():
            safe_name = sanitize_filename(trna_name)
            svg_path = org_dir / f"{safe_name}.svg"
            json_path = org_dir / f"{safe_name}.json"

            print(f"  Generating {safe_name}...")

            success = False
            for i, trna_data in enumerate(trna_copies):
                try:
                    svg_content = run_r2r(
                        trna_name,
                        trna_data["sequence"],
                        trna_data["structure"],
                    )

                    if svg_content:
                        svg_path.write_text(svg_content)

                        metadata = parse_r2r_svg(svg_path)
                        metadata["trna_name"] = trna_name
                        metadata["sequence"] = trna_data["sequence"]
                        metadata["structure"] = trna_data["structure"]
                        write_metadata(metadata, json_path)
                        n = len(metadata["nucleotides"])
                        suffix = f" (copy {i + 1})" if i > 0 else ""
                        print(f"    OK ({n} nucleotides){suffix}")
                        success = True
                        break

                except Exception as e:
                    print(f"    ERROR on copy {i + 1}: {e}")

            if not success:
                print("    FAILED: R2R could not render any copy")

    print(f"\nDone. Output in {OUTPUT_DIR}")


def extract_trnas_from_archive(archive_path: Path, ss_member: str) -> dict:
    """Extract tRNA data from a local gtRNAdb tar.gz archive.

    Returns dict mapping tRNA names to {"sequence": ..., "structure": ...}.
    """
    if not archive_path.exists():
        print(f"  Archive not found: {archive_path}")
        return {}

    print(f"  Extracting {ss_member} from {archive_path.name}")

    with tarfile.open(archive_path, "r:gz") as tar:
        member = tar.getmember(ss_member)
        f = tar.extractfile(member)
        if f is None:
            print(f"  Could not extract {ss_member}")
            return {}
        ss_data = f.read().decode("utf-8")

    return parse_trnascan_ss(ss_data)


def parse_trnascan_ss(ss_data: str) -> dict:
    """Parse tRNAscan-SE .ss output format.

    Each entry looks like:
        chr.trna71 (2518231-2518156)\tLength: 76 bp
        Type: Ala\tAnticodon: GGC at 34-36 (...)  Score: 75.0
        ...
        Seq: GGGCGTGTGGCGTAGTCGG...
        Str: >>>>>.>..>>>>........<<<<.>>>>>...

    The Str line uses > for 5' base pairs and < for 3' base pairs.
    """
    trnas = {}
    # Split on lines starting with a non-whitespace entry name
    # Handles both chr.trna* and chrN.trna* formats
    entries = re.split(r"\n(?=\S+\.trna\d+\s)", ss_data)

    for entry in entries:
        if not entry.strip():
            continue

        # Skip intron-containing tRNAs (for now)
        if "Possible intron:" in entry:
            continue

        # Extract type and anticodon
        type_match = re.search(r"Type:\s+(\S+)", entry)
        ac_match = re.search(r"Anticodon:\s+(\S+)", entry)

        if not type_match or not ac_match:
            continue

        aa_type = type_match.group(1)
        anticodon = ac_match.group(1)
        trna_name = f"tRNA-{aa_type}-{anticodon}"

        # Extract sequence and structure
        seq_match = re.search(r"Seq:\s+(\S+)", entry)
        str_match = re.search(r"Str:\s+(\S+)", entry)

        if not seq_match or not str_match:
            continue

        sequence = seq_match.group(1).upper().replace("T", "U")
        structure_raw = str_match.group(1)

        # Convert tRNAscan-SE bracket notation to R2R format
        structure = convert_trnascan_to_r2r(structure_raw)

        # Keep all occurrences; first will be tried first
        if trna_name not in trnas:
            trnas[trna_name] = []
        trnas[trna_name].append({
            "sequence": sequence,
            "structure": structure,
        })

    return trnas


def convert_trnascan_to_r2r(raw: str) -> str:
    """Convert tRNAscan-SE structure notation to R2R bracket format.

    tRNAscan uses > for 5' strand of stems and < for 3' strand.
    R2R uses < for 5' (opening) and > for 3' (closing) — the opposite.
    """
    result = []
    for ch in raw:
        if ch == ">":
            result.append("<")
        elif ch == "<":
            result.append(">")
        else:
            result.append(".")
    return "".join(result)


def extract_isotype(trna_name: str) -> str | None:
    """Extract amino acid isotype from tRNA name."""
    match = re.match(r"tRNA-(\w+)-", trna_name)
    return match.group(1) if match else None


def get_modomics_isotypes(org_name: str) -> set | None:
    """Get set of tRNA isotypes with MODOMICS modification data.

    Reads the cached MODOMICS RDS files to determine which amino acid
    types have known modifications. Returns None if no cached data.
    """
    rds_path = MODOMICS_DIR / f"{org_name}.rds"
    if not rds_path.exists():
        return None

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

    Uses the two-step R2R workflow:
    1. --GSC-weighted-consensus to build a consensus alignment
    2. --disable-usage-warning to render the SVG
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)

        # Write Stockholm file
        sto_path = tmpdir / "trna.sto"
        sto_content = create_stockholm(trna_name, sequence, structure)
        sto_path.write_text(sto_content)

        # Step 1: Build consensus
        cons_path = tmpdir / "trna.cons.sto"
        try:
            result = subprocess.run(
                [
                    R2R_BIN, "--GSC-weighted-consensus",
                    str(sto_path), str(cons_path),
                    "3", "0.97", "0.9", "0.75",
                    "4", "0.97", "0.9", "0.75", "0.5", "0.1",
                ],
                capture_output=True,
                text=True,
                timeout=60,
            )

            if result.returncode != 0:
                print(f"    R2R consensus stderr: {result.stderr[:200]}")
                return None

            if not cons_path.exists():
                print("    R2R consensus produced no output file")
                return None

        except FileNotFoundError:
            print(f"    ERROR: {R2R_BIN} not found on PATH")
            return None
        except subprocess.TimeoutExpired:
            print("    ERROR: R2R consensus timed out")
            return None

        # Step 2: Render SVG
        svg_path = tmpdir / "trna.svg"
        try:
            result = subprocess.run(
                [
                    R2R_BIN, "--disable-usage-warning",
                    str(cons_path), str(svg_path),
                ],
                capture_output=True,
                text=True,
                timeout=60,
            )

            if result.returncode != 0:
                print(f"    R2R render stderr: {result.stderr[:200]}")
                return None

            if svg_path.exists():
                return svg_path.read_text()
            return None

        except subprocess.TimeoutExpired:
            print("    ERROR: R2R render timed out")
            return None


def ensure_cca_tail(sequence: str, structure: str) -> tuple[str, str]:
    """Append CCA tail and amino acid placeholder to sequence if missing.

    Most eukaryotic tRNAs in genomic sequences lack the post-transcriptionally
    added CCA tail. We add it plus a lowercase 'a' placeholder for the amino
    acid attachment site. This extra position can be used to visualize
    aminoacylation state alongside modification data.
    """
    upper_seq = sequence.upper()
    if upper_seq.endswith("CCA"):
        # Already has CCA, just add amino acid placeholder
        sequence = sequence + "a"
        structure = structure + "."
        return sequence, structure

    # Append CCA + amino acid placeholder as unpaired nucleotides
    sequence = sequence + "CCAa"
    structure = structure + "...."
    return sequence, structure


def find_innermost_acceptor_pair(structure: str) -> int:
    """Find the position of the innermost acceptor stem base pair.

    The acceptor stem is the outermost stem starting at position 0. We walk
    inward along the 5' strand until we reach the last opening bracket
    before the multi-stem junction begins. This handles both standard
    tRNAs (7 bp acceptor) and non-standard ones (variable acceptor length).
    """
    # Build pair map
    stack = []
    pairs = {}
    for i, ch in enumerate(structure):
        if ch == "<":
            stack.append(i)
        elif ch == ">":
            if stack:
                j = stack.pop()
                pairs[j] = i
                pairs[i] = j

    # Walk inward from position 0 along the acceptor stem
    # The acceptor stem is contiguous opening brackets (possibly with
    # small internal loops of 1-2 unpaired bases)
    last_opening = None
    i = 0
    while i < len(structure):
        if structure[i] == "<":
            last_opening = i
            # Check if the region between this pair contains stems
            if i in pairs:
                partner = pairs[i]
                # Look ahead: count stems inside this pair
                inner_stems = 0
                j = i + 1
                while j < partner:
                    if structure[j] == "<" and j in pairs:
                        inner_stems += 1
                        j = pairs[j] + 1
                    else:
                        j += 1

                if inner_stems >= 3:
                    # This is the innermost acceptor pair enclosing junction
                    return i

            i += 1
        elif structure[i] == ".":
            # Allow small gaps (internal loops) in the acceptor stem
            i += 1
        else:
            break

    return last_opening if last_opening is not None else 0


def find_internal_stem_starts(structure: str, m_pos: int) -> list[int]:
    """Find the start positions of each stem in the multistem junction.

    Starting from the m_pos (innermost acceptor pair), walk the junction
    and return the position of the first opening bracket of each internal
    stem. Standard tRNAs have 3 stems (D, AC, T) while long variable arm
    tRNAs have 4 (D, AC, variable arm, T).
    """
    stack = []
    pairs = {}
    for i, ch in enumerate(structure):
        if ch == "<":
            stack.append(i)
        elif ch == ">":
            if stack:
                j = stack.pop()
                pairs[j] = i
                pairs[i] = j

    if m_pos not in pairs:
        return []

    m_pair = pairs[m_pos]
    stem_starts = []
    i = m_pos + 1
    while i < m_pair:
        if structure[i] == "<":
            stem_starts.append(i)
            if i in pairs:
                i = pairs[i] + 1
            else:
                i += 1
        else:
            i += 1

    return stem_starts


def create_stockholm(name: str, sequence: str, structure: str) -> str:
    """Create a Stockholm alignment file for R2R with layout directives.

    R2R expects standard Stockholm format with SS_cons annotation.
    Uses multistem_junction_bulgey with explicit coordinates to produce
    a classic textbook cloverleaf with straight stems:
    - Acceptor stem vertical at top
    - D-arm horizontal left
    - Anticodon arm vertical at bottom
    - T-arm horizontal right

    Labels placed:
    - ``j`` at the innermost acceptor pair (junction entry point)
    - ``3`` at the T-stem start (for place_explicit positioning)
    """
    # Ensure CCA tail and amino acid placeholder are present
    sequence, structure = ensure_cca_tail(sequence, structure)

    seq_len = len(sequence)
    str_len = len(structure)

    if str_len < seq_len:
        structure = structure + "." * (seq_len - str_len)
    elif str_len > seq_len:
        structure = structure[:seq_len]

    seq_display = sequence.replace("U", "T").replace("u", "t")
    seq_id = sanitize_filename(name)
    tag_ss = "#=GC SS_cons"
    tag_lbl = "#=GC R2R_LABEL"

    # Pad all tags to the same width
    width = max(len(seq_id), len(tag_ss), len(tag_lbl)) + 2

    # Find label positions: j at innermost acceptor pair, 1/2/3 at stems
    j_pos = find_innermost_acceptor_pair(structure)
    stem_starts = find_internal_stem_starts(structure, j_pos)
    n_stems = len(stem_starts)

    # Build label string: j at junction, numbered labels at each stem
    label_chars = ["."] * seq_len
    label_chars[j_pos] = "j"
    for idx, pos in enumerate(stem_starts):
        label_chars[pos] = str(idx + 1)
    label = "".join(label_chars)

    # R2R directives for straight-stem cloverleaf layout
    directives = []

    # Initial stem direction: 90° = vertical (acceptor points up)
    directives.append("#=GF R2R set_dir pos0 90 f")

    # Explicitly position T-stem (always the last numbered stem)
    t_label = str(n_stems)
    directives.append(
        f"#=GF R2R place_explicit {t_label} {t_label}-- 0 1 0 0 0 0"
    )

    if n_stems == 4:
        # Long variable arm (Leu/Ser): D, AC, variable arm, T
        # Angles swapped vs R2R demo so D-arm is LEFT and T-arm is RIGHT
        directives.append(
            "#=GF R2R multistem_junction_bulgey j "
            "disable_auto_flip_place_explicit "
            "J0/base 0 1.62706 2.29987 0 0 -90 "
            "J1/base 0 4.63116 0.91097 0 0 0 "
            "J2/base 0 3.5 -1.5 0 0 45 "
            "J3/base 0 2.32414 -2.32817 0 0 -270 "
            "backbonelen 1 1"
        )
    else:
        # Standard 3-stem: D, AC, T
        # Based on R2R demo coordinates with J0/J2 angles swapped
        # so D-arm points LEFT and T-arm points RIGHT
        directives.append(
            "#=GF R2R multistem_junction_bulgey j "
            "disable_auto_flip_place_explicit "
            "J0/base 0 1.62706 2.29987 0 0 -90 "
            "J1/base 0 4.63116 0.91097 0 0 0 "
            "J2/base 0 2.32414 -2.32817 0 0 -270 "
            "backbonelen 1 1"
        )

    lines = [
        "# STOCKHOLM 1.0",
        "",
        f"{seq_id:<{width}}{seq_display}",
        f"{tag_ss:<{width}}{structure}",
        f"{tag_lbl:<{width}}{label}",
        "",
        *directives,
        "",
        "//",
    ]
    return "\n".join(lines)


if __name__ == "__main__":
    main()
