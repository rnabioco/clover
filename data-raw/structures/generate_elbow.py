"""Generate tRNA elbow-form SVGs from existing cloverleaf JSON metadata.

Reads each E. coli tRNA's existing JSON (for sequence) and the Sprinzl
coordinates TSV (for per-position canonical Sprinzl labels), then emits an
elbow-form SVG and JSON in the SAME schema as the cloverleaf path. The
R-side overlay code in plot_tRNA_structure() consumes both layouts
identically.

The elbow form (per Matsumoto et al., JBC 2026; see also classical
structural-biology depictions) lays out the molecule so that the acceptor
stem and T-arm are stacked coaxially as one horizontal helix at the top,
the D-arm extends horizontally to the left below it, the anticodon arm
hangs vertically below, and the variable region forms a small bulge
between the AC-stem and T-stem.

Skips variable-arm tRNAs (Leu, Ser, SeC, Tyr, some Thr) — those will need
a v2 layout that handles the variable arm.

Usage
-----
    pixi run -e struct python data-raw/structures/generate_elbow.py
"""

from __future__ import annotations

import csv
import gzip
import json
import math
import re
import sys
from pathlib import Path

# --- paths --------------------------------------------------------------

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
SCRIPT_DIR = Path(__file__).resolve().parent
STRUCTURES_DIR = PROJECT_ROOT / "inst" / "extdata" / "structures"
SPRINZL_DIR = PROJECT_ROOT / "inst" / "extdata" / "sprinzl"

# Map organism dir name -> Sprinzl coords filename
ORG_SPRINZL = {
    "Escherichia_coli": "ecoliK12_global_coords.tsv.gz",
}

# --- layout constants ---------------------------------------------------

S = 7.56  # backbone spacing along a strand
P = 7.56  # vertical pair spacing between rows of a stem

# Combined acceptor + T helix at top
TOP_Y = 40.0          # 3' strand row (acceptor 3' + T-stem 3'); CCA emerges UP from here
BOT_Y = TOP_Y + P     # 5' strand row (acceptor 5' + T-stem 5'); side arms emerge DOWN from here
ACC_X1 = 220.0        # x of pos 1 (5' free end) and pos 72 (paired with pos 1)
GAP = S * 1.6         # break between T-stem and acceptor on the bottom strand

# T-loop (54..60) curls at far LEFT of the combined helix. With N residues
# placed on a 180° arc at t = (k+0.5)/N, neighbors are 2R*sin(π/(2N)) apart.
# Solve for R so chord = S → R = S / (2 sin(π/(2N))). Matches the cloverleaf
# spacing so mod circles (radius 4.3) fit cleanly.
T_LOOP_RADIUS = S / (2 * math.sin(math.pi / 14))  # 7 residues

# Both D-arm and AC-arm hang DOWN as parallel vertical helices below the top
# stack — D-arm to the LEFT (smaller, 4 bp + D-loop at the bottom-left),
# AC-arm to the RIGHT (larger, 5 bp + AC-loop at the bottom). Together with
# the horizontal acceptor+T helix at the top this gives the canonical L-shape
# elbow projection seen in textbook tRNA depictions (e.g., Matsumoto 2026
# JBC fig 2A).
#
# Spacing is generous enough that the radius-4.3 modification circles and
# radius ~5 outline circles used by plot_tRNA_structure() do not collide
# with neighboring nucleotide letters.
D_X_LEFT = 130.0       # left column of D-stem (pos 10-13)
D_X_RIGHT = D_X_LEFT + P
D_TOP_Y = BOT_Y + 50.0
D_LOOP_RADIUS = S / (2 * math.sin(math.pi / 16))  # 8 residues on 180° arc

AC_X_LEFT = 158.0
AC_X_RIGHT = AC_X_LEFT + P
AC_TOP_Y = D_TOP_Y     # AC-arm aligned vertically with D-arm
AC_LOOP_RADIUS = S / (2 * math.sin(math.pi / 14))  # 7 residues on 180° arc

# --- variable-arm detection --------------------------------------------

# A sprinzl label like "e1", "e12", "e23" indicates a variable-arm position;
# these tRNAs (Leu, Ser, SeC, some Tyr/Thr) have an extended variable arm
# that the v1 elbow layout does not handle. Skip them.
VAR_ARM_RE = re.compile(r"^e\d")


def main() -> int:
    total_written = 0
    total_skipped = 0

    for org_dirname, sprinzl_fname in ORG_SPRINZL.items():
        org_dir = STRUCTURES_DIR / org_dirname
        if not org_dir.is_dir():
            print(f"  SKIP {org_dirname}: structures dir not found")
            continue

        sprinzl_path = SPRINZL_DIR / sprinzl_fname
        if not sprinzl_path.exists():
            print(f"  SKIP {org_dirname}: sprinzl file not found at {sprinzl_path}")
            continue

        print(f"\n=== {org_dirname} ===")
        sprinzl = load_sprinzl(sprinzl_path)
        out_dir = org_dir / "elbow"
        out_dir.mkdir(parents=True, exist_ok=True)

        # Clean stale output
        for old in out_dir.glob("*.svg"):
            old.unlink()
        for old in out_dir.glob("*.json"):
            old.unlink()

        for json_path in sorted(org_dir.glob("*.json")):
            trna_name = json_path.stem
            try:
                metadata = json.loads(json_path.read_text())
            except Exception as e:
                print(f"  ERROR reading {trna_name}: {e}")
                continue

            sprinzl_rows = lookup_sprinzl(sprinzl, trna_name)
            if sprinzl_rows is None:
                print(f"  SKIP {trna_name}: no Sprinzl mapping found")
                total_skipped += 1
                continue

            if any(VAR_ARM_RE.match(r["sprinzl_label"]) for r in sprinzl_rows):
                print(f"  SKIP {trna_name}: variable-arm tRNA (deferred to v2)")
                total_skipped += 1
                continue

            try:
                elbow_meta = build_elbow(metadata, sprinzl_rows, trna_name)
            except KeyError as e:
                print(f"  SKIP {trna_name}: missing layout key {e}")
                total_skipped += 1
                continue

            svg_path = out_dir / f"{trna_name}.svg"
            elbow_json_path = out_dir / f"{trna_name}.json"
            svg_path.write_text(render_svg(elbow_meta, trna_name))
            elbow_json_path.write_text(json.dumps(elbow_meta, indent=2))
            n = len(elbow_meta["nucleotides"])
            print(f"  OK   {trna_name} ({n} nt)")
            total_written += 1

    print(f"\nWrote {total_written} elbow tRNA(s); skipped {total_skipped}.")
    return 0


# --- Sprinzl loading + lookup ------------------------------------------


def load_sprinzl(path: Path) -> dict[str, list[dict]]:
    """Group Sprinzl rows by trna_id."""
    by_id: dict[str, list[dict]] = {}
    with gzip.open(path, "rt") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            tid = row["trna_id"]
            try:
                row["seq_index"] = int(row["seq_index"])
            except (KeyError, ValueError):
                continue
            by_id.setdefault(tid, []).append(row)
    return by_id


def lookup_sprinzl(
    sprinzl: dict[str, list[dict]], trna_name: str
) -> list[dict] | None:
    """Find a Sprinzl row block whose trna_id matches the structure name.

    Structure names use DNA anticodon notation (e.g., 'tRNA-Asp-GTC') while
    Sprinzl files use RNA notation ('tRNA-Asp-GUC'). Convert T→U in the
    anticodon segment before matching.
    """
    rna_name = dna_to_rna_anticodon(trna_name)

    candidates = [
        tid for tid in sprinzl
        if tid == rna_name or tid.startswith(rna_name + "-")
    ]
    if not candidates:
        return None
    candidates.sort()
    rows = sprinzl[candidates[0]]
    rows.sort(key=lambda r: r["seq_index"])
    return rows


def dna_to_rna_anticodon(trna_name: str) -> str:
    m = re.match(r"^(tRNA-[A-Za-z0-9]+)-([ACGTU]+)$", trna_name)
    if not m:
        return trna_name
    aa, codon = m.group(1), m.group(2)
    return f"{aa}-{codon.replace('T', 'U')}"


# --- elbow geometry -----------------------------------------------------


def build_elbow(
    cloverleaf_meta: dict, sprinzl_rows: list[dict], trna_name: str
) -> dict:
    """Compute elbow (x, y) for each nucleotide from its Sprinzl label.

    The cloverleaf metadata gives us the per-position bases and base-pair
    line list (which we recompute against the new coords). The Sprinzl
    rows give us the Sprinzl label and region per position. We index
    nucleotides in the cloverleaf JSON by `pos` (1-based, document order)
    and join to Sprinzl rows by `seq_index`.
    """
    base_by_pos = {nuc["pos"]: nuc["base"] for nuc in cloverleaf_meta["nucleotides"]}
    sprinzl_by_pos = {r["seq_index"]: r for r in sprinzl_rows}

    coord_table = elbow_coord_table()

    nucleotides = []
    pos_to_xy: dict[int, tuple[float, float]] = {}
    for pos in sorted(base_by_pos.keys()):
        spr = sprinzl_by_pos.get(pos)
        if spr is None:
            raise KeyError(f"pos {pos} missing from Sprinzl rows")
        label = spr["sprinzl_label"]
        if label not in coord_table:
            raise KeyError(f"label {label!r} (pos {pos})")
        x, y = coord_table[label]
        pos_to_xy[pos] = (x, y)
        nucleotides.append({
            "base": base_by_pos[pos],
            "x": round(x, 2),
            "y": round(y, 2),
            "pos": pos,
        })

    lines = build_pair_lines(sprinzl_by_pos, pos_to_xy)

    xs = [n["x"] for n in nucleotides]
    ys = [n["y"] for n in nucleotides]
    width = round(max(xs) + 20, 2)
    height = round(max(ys) + 20, 2)

    return {
        "nucleotides": nucleotides,
        "width": width,
        "height": height,
        "lines": lines,
        "trna_name": trna_name,
        "sequence": cloverleaf_meta.get("sequence", ""),
        "structure": cloverleaf_meta.get("structure", ""),
    }


def elbow_coord_table() -> dict[str, tuple[float, float]]:
    """Sprinzl label -> (x, y) for the elbow layout.

    Top combined helix:
      Top row (3' strand, continuous backbone L->R):
        61 62 63 64 65 66 67 68 69 70 71 72
      Bottom row (5' strand, with gap between T-stem and acceptor):
        53 52 51 50 49  [GAP]  7 6 5 4 3 2 1
      Pair vertically: 1<->72, 2<->71, ... 7<->66; 49<->65, 50<->64, ... 53<->61.

    T-loop (54..60) curls at far LEFT of the combined helix.
    D-arm sits BELOW the top stack, horizontal, with D-loop at far LEFT.
    AC-arm hangs DOWN below the rest, with AC-loop at the BOTTOM.
    Variable region (44..48) curves from the right of the AC-stem 3' end
    UP and OVER to position 49 on the top stack.
    Discriminator + CCA (73..76) extend UP from pos 72 (top-right).
    """
    coords: dict[str, tuple[float, float]] = {}

    # ----- acceptor stem (1..7 bottom row, 66..72 top row) -----
    for n in range(1, 8):
        x = ACC_X1 - (n - 1) * S
        coords[str(n)] = (x, BOT_Y)
    for n in range(66, 73):
        x = ACC_X1 - (72 - n) * S
        coords[str(n)] = (x, TOP_Y)

    # rightmost / leftmost x of acceptor section
    acc_left_x = ACC_X1 - 6 * S        # x of pos 7 / pos 66
    t_right_x = acc_left_x - GAP        # x of pos 49 / pos 65 (across the gap)

    # ----- T-stem (49..53 bottom row, 61..65 top row) -----
    for n in range(49, 54):
        x = t_right_x - (n - 49) * S
        coords[str(n)] = (x, BOT_Y)
    for n in range(61, 66):
        x = t_right_x - (65 - n) * S
        coords[str(n)] = (x, TOP_Y)

    t_left_x = t_right_x - 4 * S  # x of pos 53 / pos 61

    # ----- T-loop (54..60), 7 residues on a 180° arc curling LEFT -----
    arc_cx = t_left_x - T_LOOP_RADIUS * 0.35
    arc_cy = (TOP_Y + BOT_Y) / 2
    for k, n in enumerate(range(54, 61)):
        # k = 0..6, sweep from "just past pos 53" (bottom) up & around to
        # "just before pos 61" (top). In SVG y-down: angle pi/2 = down,
        # angle 3pi/2 = up; sweeping through pi (left).
        t = (k + 0.5) / 7  # 1/14, 3/14, ..., 13/14
        theta = math.pi / 2 + t * math.pi
        x = arc_cx + T_LOOP_RADIUS * math.cos(theta)
        y = arc_cy + T_LOOP_RADIUS * math.sin(theta)
        coords[str(n)] = (x, y)

    # ----- discriminator + CCA (73..76) extend UP from pos 72 -----
    cca_x = ACC_X1
    for k, n in enumerate(range(73, 77)):
        # 73 just above 72 on the same column, 74-76 stacking up
        coords[str(n)] = (cca_x, TOP_Y - (k + 1) * S)

    # ----- D-arm: VERTICAL stem, D-loop curls at BOTTOM-LEFT -----
    # pos 10 (D-stem 5' first) at TOP-LEFT — paired with pos 25 at TOP-RIGHT.
    # Stem extends DOWN from there; pos 13 at BOTTOM-LEFT, pos 22 at BOTTOM-
    # RIGHT. D-loop (14..21) is a 180° arc connecting pos 13 (bottom-left)
    # around the LEFT side back up to pos 22 (bottom-right) — i.e., it
    # bulges LEFT and DOWN from the bottom of the stem, matching the
    # canonical L-shape projection.
    # Hinges 8, 9 are the diagonal connector from pos 7 (top stack, bot row,
    # at acc_left_x) DOWN-LEFT to pos 10 (top-left of D-arm). Place evenly
    # along the line so consecutive backbone neighbors are similarly spaced.
    for k, n in enumerate(("8", "9"), start=1):
        t = k / 3.0  # 1/3 and 2/3 along the line from pos 7 to pos 10
        coords[n] = (
            acc_left_x + (D_X_LEFT - acc_left_x) * t,
            BOT_Y + (D_TOP_Y - BOT_Y) * t,
        )
    # D-stem 5' (10..13) on LEFT column, top-to-bottom
    for n in range(10, 14):
        coords[str(n)] = (D_X_LEFT, D_TOP_Y + (n - 10) * S)
    # D-stem 3' (22..25) on RIGHT column; pair: 25-10, 24-11, 23-12, 22-13
    for n in range(22, 26):
        coords[str(n)] = (D_X_RIGHT, D_TOP_Y + (25 - n) * S)

    d_bot_y = D_TOP_Y + 3 * S  # y of pos 13 / pos 22

    # ----- D-loop (14..21), 8 residues on a 180° arc curling LEFT-and-DOWN -----
    # Arc center sits a bit DOWN-LEFT of the stem bottom. Sweep angle covers
    # the LEFT half of the circle, from "just past pos 13" (top of arc, bot
    # of stem on left col) curving DOWN through far-LEFT and back UP to
    # "just before pos 22" (top of arc, bot of stem on right col).
    arc_cx_d = D_X_LEFT + P / 2
    arc_cy_d = d_bot_y + D_LOOP_RADIUS * 0.55
    for k, n in enumerate(range(14, 22)):
        t = (k + 0.5) / 8
        # Parameterize: angle from pi (left of stem-left-column-bottom)
        # sweeping CCW through 3pi/2 (BOTTOM of arc) to 2pi (right side).
        # In SVG y-down, y = cy + r*sin(theta) — sin pi to 2pi gives 0 to 0
        # via -1 at 3pi/2; we want POSITIVE y at the bottom of the arc, so
        # use abs offset.
        theta = math.pi + t * math.pi
        x = arc_cx_d + D_LOOP_RADIUS * math.cos(theta)
        y = arc_cy_d - D_LOOP_RADIUS * math.sin(theta)
        coords[str(n)] = (x, y)

    # D-loop insertions: 17a is between 17 and 18; 20a between 20 and 21;
    # 20b between 20a and 21. Place them on the same arc with offset t.
    insertion_ts = {
        "17a": (3.5 + 0.5) / 8,
        "20a": (6.3 + 0.5) / 8,
        "20b": (6.6 + 0.5) / 8,
    }
    for label, t in insertion_ts.items():
        theta = math.pi + t * math.pi
        x = arc_cx_d + D_LOOP_RADIUS * math.cos(theta)
        y = arc_cy_d - D_LOOP_RADIUS * math.sin(theta)
        coords[label] = (x, y)

    # ----- pos 26 (hinge between D-stem 3' top and AC-stem 5' top) -----
    # short horizontal bridge from D-arm TOP-RIGHT (pos 25) over to AC-arm
    # TOP-LEFT (pos 27). Both at the same y as the top of the side arms.
    coords["26"] = (
        (D_X_RIGHT + AC_X_LEFT) / 2,
        D_TOP_Y - 2.5,
    )

    # ----- anticodon stem (27..31 left col, 39..43 right col) -----
    # VERTICAL helix to the RIGHT of the D-arm. Pairs: 27-43, 28-42, ..., 31-39
    # 27 at top-left, 31 at bottom-left (just above AC-loop).
    for n in range(27, 32):
        y = AC_TOP_Y + (n - 27) * S
        coords[str(n)] = (AC_X_LEFT, y)
    for n in range(39, 44):
        y = AC_TOP_Y + (43 - n) * S
        coords[str(n)] = (AC_X_RIGHT, y)

    ac_bot_y = AC_TOP_Y + 4 * S  # y of pos 31 / pos 39

    # ----- AC-loop (32..38), 7 residues on 180° arc curling DOWN -----
    arc_cx_ac = (AC_X_LEFT + AC_X_RIGHT) / 2
    arc_cy_ac = ac_bot_y + AC_LOOP_RADIUS * 0.35
    for k, n in enumerate(range(32, 39)):
        t = (k + 0.5) / 7
        # sweep from angle pi (left, near pos 31 column) to 0 (right, near pos 39 column)
        # going DOWN through 3pi/2 (in y-down convention, that's bottom).
        # parameterize: theta = pi - t * pi  (pi to 0)
        # but we want the arc to bulge DOWNWARD (positive y in SVG)
        theta = math.pi - t * math.pi
        x = arc_cx_ac + AC_LOOP_RADIUS * math.cos(theta)
        # bulge down: use abs(sin) so all loop residues have y > arc_cy_ac
        y = arc_cy_ac + AC_LOOP_RADIUS * abs(math.sin(theta))
        coords[str(n)] = (x, y)

    # ----- variable region (44..48): smooth curve from AC-stem 3' top to pos 49 -----
    # The bulge runs from pos 43 (top-right of AC-stem) up to pos 49 (T-stem
    # 5' end on bottom row of top stack). Quadratic bezier with control
    # point bowed RIGHT so the path doesn't crowd the AC-stem 3' column.
    # Positions are placed by ARC-LENGTH along the curve so consecutive
    # neighbors are evenly spaced (a t-uniform sampling clusters at the apex).
    var_start = (AC_X_RIGHT, AC_TOP_Y - 0.6 * S)
    var_end = (t_right_x + 0.3 * S, BOT_Y + 0.6 * S)
    ctrl_x = max(var_start[0], var_end[0]) + 18
    ctrl_y = (var_start[1] + var_end[1]) / 2

    def bez(t: float) -> tuple[float, float]:
        x = (1 - t) ** 2 * var_start[0] + 2 * (1 - t) * t * ctrl_x + t ** 2 * var_end[0]
        y = (1 - t) ** 2 * var_start[1] + 2 * (1 - t) * t * ctrl_y + t ** 2 * var_end[1]
        return x, y

    # Sample the curve densely to build an arc-length lookup
    n_samples = 200
    samples = [bez(i / n_samples) for i in range(n_samples + 1)]
    cum_lengths = [0.0]
    for a, b in zip(samples, samples[1:]):
        cum_lengths.append(cum_lengths[-1] + math.hypot(b[0] - a[0], b[1] - a[1]))
    total_length = cum_lengths[-1]

    def t_for_arc(target_arc: float) -> float:
        for i, L in enumerate(cum_lengths):
            if L >= target_arc:
                # linear interpolate t around index i
                if i == 0:
                    return 0.0
                prev = cum_lengths[i - 1]
                frac = (target_arc - prev) / (L - prev) if L > prev else 0.0
                return ((i - 1) + frac) / n_samples
        return 1.0

    # 5 residues, evenly spaced along the curve (1/6, 2/6, ..., 5/6 of total)
    for k, n in enumerate(range(44, 49), start=1):
        target = total_length * k / 6
        t = t_for_arc(target)
        coords[str(n)] = bez(t)

    return coords


def build_pair_lines(
    sprinzl_by_pos: dict[int, dict],
    pos_to_xy: dict[int, tuple[float, float]],
) -> list[dict]:
    """Compute base-pair tick lines from Sprinzl pairings."""
    label_to_pos: dict[str, int] = {}
    for pos, row in sprinzl_by_pos.items():
        label = row["sprinzl_label"]
        label_to_pos.setdefault(label, pos)

    canonical_pairs = [
        # acceptor stem
        ("1", "72"), ("2", "71"), ("3", "70"), ("4", "69"),
        ("5", "68"), ("6", "67"), ("7", "66"),
        # D-stem
        ("10", "25"), ("11", "24"), ("12", "23"), ("13", "22"),
        # AC-stem
        ("27", "43"), ("28", "42"), ("29", "41"),
        ("30", "40"), ("31", "39"),
        # T-stem
        ("49", "65"), ("50", "64"), ("51", "63"),
        ("52", "62"), ("53", "61"),
    ]

    lines: list[dict] = []
    for a, b in canonical_pairs:
        pa, pb = label_to_pos.get(a), label_to_pos.get(b)
        if pa is None or pb is None:
            continue
        xa, ya = pos_to_xy[pa]
        xb, yb = pos_to_xy[pb]
        # base-pair tick: short segment between the two paired letter centers,
        # offset so it doesn't overlap the letters.
        dx, dy = xb - xa, yb - ya
        L = math.hypot(dx, dy)
        if L < 1e-6:
            continue
        ux, uy = dx / L, dy / L
        pad = 2.5
        x1 = xa + ux * pad
        y1 = ya + uy * pad
        x2 = xb - ux * pad
        y2 = yb - uy * pad
        # add small character-center offset (matches plot-structure.R nuc_*_offset)
        x1 += 2.375
        y1 += -2.565
        x2 += 2.375
        y2 += -2.565
        lines.append({
            "x1": round(x1, 2),
            "y1": round(y1, 2),
            "x2": round(x2, 2),
            "y2": round(y2, 2),
        })
    return lines


# --- SVG rendering ------------------------------------------------------


SVG_HEADER = """<svg
   xmlns:svg="http://www.w3.org/2000/svg"
   xmlns="http://www.w3.org/2000/svg"
   xmlns:xlink="http://www.w3.org/1999/xlink"
   version="1.0"
   width="{w}"
   height="{h}"
   xml:space="preserve">
<text x="0" y="8.616">
  <tspan x="0" y="8.616" fill="#000000" font-variant="normal" font-weight="normal" font-style="normal" font-family="Helvetica, Arial, sans-serif" font-size="9">{name}</tspan>
</text>
"""

NUC_TEXT = (
    '<text x="{x}" y="{y}">\n'
    '  <tspan x="{x}" y="{y}" fill="#d90000" font-variant="normal" '
    'font-weight="normal" font-style="normal" '
    'font-family="Helvetica, Arial, sans-serif" font-size="7.1">{base}</tspan>\n'
    '</text>\n'
)

PAIR_PATH = (
    '<path fill="none" stroke="#000000" stroke-width="1.44" '
    'stroke-linecap="butt" stroke-linejoin="miter" '
    'd="M {x1} {y1} L {x2} {y2}"/>\n'
)


def render_svg(meta: dict, trna_name: str) -> str:
    parts = [SVG_HEADER.format(w=meta["width"], h=meta["height"], name=trna_name)]
    for line in meta["lines"]:
        parts.append(PAIR_PATH.format(**{
            "x1": line["x1"], "y1": line["y1"],
            "x2": line["x2"], "y2": line["y2"],
        }))
    for nuc in meta["nucleotides"]:
        parts.append(NUC_TEXT.format(x=nuc["x"], y=nuc["y"], base=nuc["base"]))
    parts.append("</svg>\n")
    return "".join(parts)


if __name__ == "__main__":
    sys.exit(main())
