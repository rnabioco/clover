"""Parse R2R SVG output to extract nucleotide positions and coordinates.

Produces a JSON metadata file mapping each nucleotide to its SVG (x, y)
coordinates, sequence position, and base identity. This metadata is used
by the R runtime functions to overlay modification annotations.
"""

import json
import re
from pathlib import Path

from lxml import etree

SVG_NS = "http://www.w3.org/2000/svg"
NSMAP = {"svg": SVG_NS}


def parse_r2r_svg(svg_path: str | Path) -> dict:
    """Parse an R2R-generated SVG and extract nucleotide metadata.

    R2R renders each nucleotide as a <text> element with a single
    character (A, C, G, U, or T) at a specific (x, y) position.

    Parameters
    ----------
    svg_path : str or Path
        Path to the R2R SVG file.

    Returns
    -------
    dict with keys:
        - "nucleotides": list of dicts with keys "pos", "base", "x", "y"
        - "width": SVG width
        - "height": SVG height
        - "lines": list of dicts with keys "x1", "y1", "x2", "y2"
          (base-pair lines)
    """
    svg_path = Path(svg_path)
    tree = etree.parse(str(svg_path))
    root = tree.getroot()

    # Get SVG dimensions
    width = float(root.get("width", "0").replace("px", ""))
    height = float(root.get("height", "0").replace("px", ""))

    # If width/height are 0, try viewBox
    if width == 0 or height == 0:
        viewbox = root.get("viewBox", "")
        if viewbox:
            parts = viewbox.split()
            if len(parts) == 4:
                width = float(parts[2])
                height = float(parts[3])

    nucleotides = []
    lines = []

    # Find all text elements - R2R puts single nucleotide chars in <tspan>
    # elements inside <text>, not as direct text content.
    for text_el in root.iter(f"{{{SVG_NS}}}text"):
        # Check direct text content first
        content = (text_el.text or "").strip()

        # Also check child <tspan> elements (R2R's actual format)
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

    # Find all line elements (base pairs)
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

    # Sort nucleotides by reading order: top-to-bottom, left-to-right
    # This usually corresponds to 5'->3' in a cloverleaf layout
    # We assign positions after sorting
    nucleotides.sort(key=lambda n: (round(n["y"], 0), n["x"]))

    # Assign sequence positions (1-based)
    for i, nuc in enumerate(nucleotides):
        nuc["pos"] = i + 1

    return {
        "nucleotides": nucleotides,
        "width": round(width, 2),
        "height": round(height, 2),
        "lines": lines,
    }


def assign_sprinzl_labels(metadata: dict, sequence: str, structure: str) -> dict:
    """Assign Sprinzl position labels to nucleotides.

    Uses the known tRNA secondary structure to assign canonical Sprinzl
    numbering. The structure string uses dot-bracket notation where
    paired positions get stem numbers and unpaired positions get loop
    numbers.

    Parameters
    ----------
    metadata : dict
        Output from parse_r2r_svg().
    sequence : str
        tRNA sequence (should match nucleotides in SVG).
    structure : str
        Dot-bracket secondary structure string.

    Returns
    -------
    dict with updated nucleotides containing "sprinzl_pos" field.
    """
    # Standard tRNA Sprinzl positions for canonical tRNA (76 positions)
    # Acceptor stem: 1-7, 66-72
    # D-stem: 10-13, 22-25
    # D-loop: 14-21 (with possible insertions 17a, 20a, 20b)
    # Anticodon stem: 27-31, 39-43
    # Anticodon loop: 32-38
    # Variable loop: 44-48 (with possible insertions e11-e15)
    # T-stem: 49-53, 61-65
    # T-loop: 54-60
    # Discriminator: 73
    # CCA: 74-76

    # For now, assign simple 1-based positions matching the sequence
    nucs = metadata["nucleotides"]
    for i, nuc in enumerate(nucs):
        if i < len(sequence):
            nuc["sprinzl_pos"] = str(i + 1)
        else:
            nuc["sprinzl_pos"] = ""

    return metadata


def write_metadata(metadata: dict, output_path: str | Path) -> None:
    """Write nucleotide metadata to a JSON file."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as f:
        json.dump(metadata, f, indent=2)


def _get_x(element) -> float | None:
    """Get x coordinate from a text element, checking attributes."""
    x = element.get("x")
    if x is not None:
        return _safe_float(x)
    # Check transform for translate
    transform = element.get("transform", "")
    match = re.search(r"translate\(\s*([\d.e+-]+)", transform)
    if match:
        return float(match.group(1))
    return None


def _get_y(element) -> float | None:
    """Get y coordinate from a text element, checking attributes."""
    y = element.get("y")
    if y is not None:
        return _safe_float(y)
    transform = element.get("transform", "")
    match = re.search(r"translate\(\s*[\d.e+-]+\s*[,\s]\s*([\d.e+-]+)", transform)
    if match:
        return float(match.group(1))
    return None


def _safe_float(val: str | None) -> float | None:
    """Safely convert a string to float."""
    if val is None:
        return None
    try:
        return float(val.replace("px", ""))
    except (ValueError, TypeError):
        return None
