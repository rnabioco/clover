#!/usr/bin/env python3
"""Re-parse all existing SVGs to regenerate JSON metadata files.

Uses the corrected document-order position assignment (no spatial sort).
Only requires Python stdlib (xml.etree.ElementTree), no lxml needed.
"""

import json
import re
from pathlib import Path
from xml.etree import ElementTree as ET

SVG_NS = "http://www.w3.org/2000/svg"
STRUCTURES_DIR = Path(__file__).resolve().parent.parent.parent / "inst" / "extdata" / "structures"


def parse_svg_stdlib(svg_path: Path) -> dict:
    """Parse an R2R SVG using stdlib ElementTree.

    Returns dict with nucleotides (in document order), width, height, lines.
    """
    tree = ET.parse(str(svg_path))
    root = tree.getroot()

    # Get SVG dimensions
    width_str = re.sub(r"[a-zA-Z]+", "", root.get("width", "0"))
    height_str = re.sub(r"[a-zA-Z]+", "", root.get("height", "0"))
    width = float(width_str) if width_str else 0
    height = float(height_str) if height_str else 0

    if width == 0 or height == 0:
        viewbox = root.get("viewBox", "")
        if viewbox:
            parts = viewbox.split()
            if len(parts) == 4:
                width = float(parts[2])
                height = float(parts[3])

    nucleotides = []
    lines = []

    for text_el in root.iter(f"{{{SVG_NS}}}text"):
        content = (text_el.text or "").strip()

        if not (len(content) == 1 and content in "ACGUTacgut"):
            for tspan in text_el.findall(f"{{{SVG_NS}}}tspan"):
                tspan_text = (tspan.text or "").strip()
                if len(tspan_text) == 1 and tspan_text in "ACGUTacgut":
                    content = tspan_text
                    break

        if len(content) == 1 and content in "ACGUTacgut":
            x = _get_x(text_el)
            y = _get_y(text_el)
            if x is not None and y is not None:
                nucleotides.append({
                    "base": content.upper(),
                    "x": round(x, 2),
                    "y": round(y, 2),
                })

    for line_el in root.iter(f"{{{SVG_NS}}}line"):
        x1 = _safe_float(line_el.get("x1"))
        y1 = _safe_float(line_el.get("y1"))
        x2 = _safe_float(line_el.get("x2"))
        y2 = _safe_float(line_el.get("y2"))
        if all(v is not None for v in [x1, y1, x2, y2]):
            lines.append({
                "x1": round(x1, 2),
                "y1": round(y1, 2),
                "x2": round(x2, 2),
                "y2": round(y2, 2),
            })

    # Document order IS the correct 5'->3' sequence order — no sorting.
    for i, nuc in enumerate(nucleotides):
        nuc["pos"] = i + 1

    return {
        "nucleotides": nucleotides,
        "width": round(width, 2),
        "height": round(height, 2),
        "lines": lines,
    }


def _get_x(element) -> float | None:
    x = element.get("x")
    if x is not None:
        return _safe_float(x)
    transform = element.get("transform", "")
    match = re.search(r"translate\(\s*([\d.e+-]+)", transform)
    if match:
        return float(match.group(1))
    return None


def _get_y(element) -> float | None:
    y = element.get("y")
    if y is not None:
        return _safe_float(y)
    transform = element.get("transform", "")
    match = re.search(r"translate\(\s*[\d.e+-]+\s*[,\s]\s*([\d.e+-]+)", transform)
    if match:
        return float(match.group(1))
    return None


def _safe_float(val: str | None) -> float | None:
    if val is None:
        return None
    try:
        return float(val.replace("px", ""))
    except (ValueError, TypeError):
        return None


def main():
    if not STRUCTURES_DIR.exists():
        print(f"ERROR: structures directory not found: {STRUCTURES_DIR}")
        return

    json_files = sorted(STRUCTURES_DIR.glob("**/tRNA-*.json"))
    print(f"Found {len(json_files)} JSON files to regenerate")

    for json_path in json_files:
        svg_path = json_path.with_suffix(".svg")
        if not svg_path.exists():
            print(f"  SKIP {json_path.name} (no matching SVG)")
            continue

        # Read existing JSON to preserve metadata fields
        with open(json_path) as f:
            old_data = json.load(f)

        # Re-parse SVG with corrected document-order logic
        new_data = parse_svg_stdlib(svg_path)

        seq = old_data.get("sequence", "")
        seq_len = len(seq)

        # Trim to sequence length (exclude amino acid placeholder)
        if seq_len > 0 and len(new_data["nucleotides"]) > seq_len:
            new_data["nucleotides"] = new_data["nucleotides"][:seq_len]

        # Verify bases match sequence
        parsed_bases = "".join(n["base"] for n in new_data["nucleotides"])
        if seq and parsed_bases != seq:
            print(f"  WARN {json_path.name}: bases mismatch!")
            print(f"    seq:    {seq}")
            print(f"    parsed: {parsed_bases}")

        # Preserve metadata from old JSON
        new_data["trna_name"] = old_data.get("trna_name", "")
        new_data["sequence"] = seq
        new_data["structure"] = old_data.get("structure", "")

        with open(json_path, "w") as f:
            json.dump(new_data, f, indent=2)
            f.write("\n")

        print(f"  OK {json_path.relative_to(STRUCTURES_DIR)}")

    print("Done.")


if __name__ == "__main__":
    main()
