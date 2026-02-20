"""Generate tRNA cloverleaf SVGs using R2R.

Aligns mature tRNA FASTA sequences against tRNAscan-SE covariance
models using Infernal's cmalign, extracts per-sequence secondary
structures, converts to Stockholm format, runs R2R to produce
cloverleaf SVGs, and extracts position metadata as JSON.

Usage:
    pixi run -e struct python data-raw/structures/generate_svgs.py

Prerequisites:
    - Python with lxml
    - R2R and cmalign on PATH (pixi struct environment)
    - Mature tRNA FASTAs in data-raw/structures/fasta/
    - tRNAscan-SE CMs in data-raw/structures/
"""

import re
import subprocess
import tempfile
from pathlib import Path

from parse_r2r_svg import parse_r2r_svg, write_metadata

# Project root (two levels up from this script)
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
SCRIPT_DIR = Path(__file__).resolve().parent
FASTA_DIR = SCRIPT_DIR / "fasta"
CM_DIR = SCRIPT_DIR
OUTPUT_DIR = PROJECT_ROOT / "inst" / "extdata" / "structures"
MODOMICS_DIR = PROJECT_ROOT / "inst" / "extdata" / "modomics"

# r2r binary — multistem_junction_bulgey doesn't need the NLOPT solver,
# so we can use whichever r2r is on PATH
R2R_BIN = "r2r"

ORGANISMS = {
    "Escherichia_coli": {
        "fasta": FASTA_DIR / "eschColi_K_12_MG1655-mature-tRNAs.fa",
        "cm": CM_DIR / "TRNAinf-bact.cm",
    },
    "Saccharomyces_cerevisiae": {
        "fasta": FASTA_DIR / "sacCer3-mature-tRNAs.fa",
        "cm": CM_DIR / "TRNAinf-euk.cm",
    },
    "Homo_sapiens": {
        "fasta": FASTA_DIR / "hg38-mature-tRNAs.fa",
        "cm": CM_DIR / "TRNAinf-euk.cm",
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

        # Step 1: Align sequences against tRNAscan-SE CM
        fasta_path = org_info["fasta"]
        if not fasta_path.exists():
            print(f"  WARNING: FASTA not found: {fasta_path}")
            continue

        cm_path = org_info["cm"]
        sto_path = run_cmalign(fasta_path, cm_path)
        if sto_path is None:
            print(f"  WARNING: cmalign failed for {org_name}")
            continue

        # Step 2: Parse aligned sequences and structures
        trnas = parse_cmalign_stockholm(sto_path)

        if not trnas:
            print(f"  WARNING: No tRNA data extracted for {org_name}")
            continue

        print(f"  Found {len(trnas)} unique tRNA isotype-anticodon combos")

        # Step 3: Clean output directory of stale files
        for old_file in org_dir.glob("*.svg"):
            old_file.unlink()
        for old_file in org_dir.glob("*.json"):
            old_file.unlink()

        # Step 4: Filter to isotypes with MODOMICS data
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

        # Step 5: Generate SVGs for each tRNA
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


def run_cmalign(fasta_path: Path, cm_path: Path) -> Path | None:
    """Run Infernal cmalign on a FASTA file against a covariance model.

    Returns path to the Stockholm output file, or None on failure.
    """
    print(f"  Running cmalign on {fasta_path.name}...")

    sto_fd, sto_path = tempfile.mkstemp(suffix=".sto")
    sto_path = Path(sto_path)

    try:
        result = subprocess.run(
            [
                "cmalign", "--notrunc",
                "-o", str(sto_path),
                str(cm_path), str(fasta_path),
            ],
            capture_output=True,
            text=True,
            timeout=300,
        )

        if result.returncode != 0:
            print(f"  cmalign stderr: {result.stderr[:500]}")
            sto_path.unlink(missing_ok=True)
            return None

        print(f"  cmalign completed successfully")
        return sto_path

    except FileNotFoundError:
        print("  ERROR: cmalign not found on PATH")
        sto_path.unlink(missing_ok=True)
        return None
    except subprocess.TimeoutExpired:
        print("  ERROR: cmalign timed out")
        sto_path.unlink(missing_ok=True)
        return None


def parse_cmalign_stockholm(sto_path: Path) -> dict:
    """Parse cmalign Stockholm output into per-tRNA sequences and structures.

    Reads aligned sequences and the SS_cons line, then for each sequence:
    1. Extract tRNA name from FASTA header (tRNA-{AA}-{Anticodon})
    2. Remove gap columns (where seq has '-' or '.')
    3. Convert WUSS notation to R2R bracket format
    4. Group by tRNA-{AA}-{Anticodon} (dedup copies)

    Returns dict mapping tRNA names to list of {sequence, structure}.
    """
    text = sto_path.read_text()

    # Parse interleaved Stockholm: collect sequences and SS_cons
    seq_parts = {}  # name -> list of seq chunks (in order)
    ss_parts = []   # list of SS_cons chunks (in order)

    for line in text.splitlines():
        line = line.rstrip()
        if not line or line.startswith("#=GF") or line.startswith("//"):
            continue
        if line.startswith("# STOCKHOLM"):
            continue
        if line.startswith("#=GC SS_cons"):
            # Extract structure after the tag
            ss_chunk = line.split(None, 2)[2] if len(line.split(None, 2)) > 2 else ""
            ss_parts.append(ss_chunk)
        elif line.startswith("#=GC") or line.startswith("#=GR"):
            continue
        elif not line.startswith("#"):
            # Sequence line: "name  aligned_seq"
            parts = line.split(None, 1)
            if len(parts) == 2:
                name, seq_chunk = parts
                if name not in seq_parts:
                    seq_parts[name] = []
                seq_parts[name].append(seq_chunk)

    if not ss_parts:
        print("  WARNING: No SS_cons found in Stockholm output")
        return {}

    # Concatenate interleaved chunks
    full_ss = "".join(ss_parts)
    full_seqs = {name: "".join(chunks) for name, chunks in seq_parts.items()}

    # Extract per-sequence ungapped structures
    trnas = {}
    name_pattern = re.compile(r"tRNA-(\w+)-(\w+)-\d+-\d+")

    for seq_name, aligned_seq in full_seqs.items():
        # Extract tRNA identity from name
        match = name_pattern.search(seq_name)
        if not match:
            print(f"  WARNING: Could not parse tRNA name from: {seq_name}")
            continue

        aa_type = match.group(1)
        anticodon = match.group(2)
        trna_name = f"tRNA-{aa_type}-{anticodon}"

        # Remove gap columns for this sequence
        ungapped_seq = []
        ungapped_ss = []
        for seq_ch, ss_ch in zip(aligned_seq, full_ss):
            if seq_ch not in ("-", "."):
                ungapped_seq.append(seq_ch)
                ungapped_ss.append(ss_ch)

        sequence = "".join(ungapped_seq).upper()
        structure = wuss_to_r2r("".join(ungapped_ss))

        if trna_name not in trnas:
            trnas[trna_name] = []
        trnas[trna_name].append({
            "sequence": sequence,
            "structure": structure,
        })

    return trnas


def wuss_to_r2r(wuss: str) -> str:
    """Convert WUSS notation (from cmalign SS_cons) to R2R bracket format.

    WUSS paired characters → R2R < and >:
        ( → <    ) → >
        < → <    > → >
        [ → <    ] → >
        { → <    } → >
    WUSS unpaired characters → R2R dot:
        : , _ - ~ . → .
    """
    opening = set("(<[{")
    closing = set(")>]}")
    unpaired = set(":,_-~.")

    result = []
    for ch in wuss:
        if ch in opening:
            result.append("<")
        elif ch in closing:
            result.append(">")
        elif ch in unpaired:
            result.append(".")
        else:
            # Unknown character (e.g., pseudoknot markers A-a, B-b)
            # Treat as unpaired for R2R
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
                return postprocess_svg(svg_path.read_text(), trna_name)
            return None

        except subprocess.TimeoutExpired:
            print("    ERROR: R2R render timed out")
            return None


def postprocess_svg(svg_content: str, trna_name: str) -> str:
    """Post-process R2R SVG output to fix label, font, and remove clutter.

    - Replace default "trna.cons" label with the actual tRNA name
    - Replace Bitstream Vera Sans font with standard sans-serif stack
    - Remove pink base-pair conservation boxes (not useful for single seqs)
    - Remove empty backbone placeholder paths
    """
    svg_content = svg_content.replace(
        ">trna.cons<", f">{trna_name}<"
    )
    svg_content = svg_content.replace(
        'font-family="Bitstream Vera Sans"',
        'font-family="Helvetica, Arial, sans-serif"',
    )

    # Remove pink base-pair boxes (fill="#ffd8d8")
    svg_content = re.sub(
        r'<path\s*\n\s*fill="#ffd8d8"[^/]*/>\n',
        "",
        svg_content,
    )

    # Remove empty backbone placeholder paths (stroke="#5c5c5c" with d="")
    svg_content = re.sub(
        r'<path\s*\n\s*fill="none"\s+stroke="#5c5c5c"[^/]*\n\s*d=""\s*/>\n',
        "",
        svg_content,
    )

    return svg_content


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
    - ``1``, ``2``, etc. at each internal stem start
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

    if n_stems == 4:
        # Long variable arm (Leu/Ser): D, AC, variable arm, T
        # Angles from R2R demo: D-arm LEFT, T-arm RIGHT
        directives.append(
            "#=GF R2R multistem_junction_bulgey j "
            "disable_auto_flip_place_explicit "
            "J0/base 0 1.62706 2.29987 0 0 -270 "
            "J1/base 0 4.63116 0.91097 0 0 0 "
            "J2/base 0 3.5 -1.5 0 0 -45 "
            "J3/base 0 2.32414 -2.32817 0 0 -90 "
            "backbonelen 1 1"
        )
    else:
        # Standard 3-stem: D, AC, T
        # Angles from R2R demo: D-arm LEFT, T-arm RIGHT
        directives.append(
            "#=GF R2R multistem_junction_bulgey j "
            "disable_auto_flip_place_explicit "
            "J0/base 0 1.62706 2.29987 0 0 -270 "
            "J1/base 0 4.63116 0.91097 0 0 0 "
            "J2/base 0 2.32414 -2.32817 0 0 -90 "
            "backbonelen 1 1"
        )

    # Set direction for discriminator base (first unpaired 3' position)
    # so the CCA tail points straight up, using turn_ss
    disc_pos = structure.rindex(">") + 1
    label_chars[disc_pos] = "d"
    label = "".join(label_chars)
    directives.append("#=GF R2R turn_ss d 0")

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
