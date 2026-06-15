"""
TWIST UDI plate well → Dorado barcode label mapping.

Layout matches Dorado TWIST-96{A,B,C,D}-UDI kits (column-major: A01, B01, … H01, A02, … H12).
"""

from __future__ import annotations

import re
from pathlib import Path

TWIST_KITS = frozenset(
    {
        "TWIST-96A-UDI",
        "TWIST-96B-UDI",
        "TWIST-96C-UDI",
        "TWIST-96D-UDI",
    }
)

PLATE_ROWS = "ABCDEFGH"
WELL_ID_RE = re.compile(r"^([A-H])(?:0?([1-9]|1[0-2]))$", re.IGNORECASE)
DORADO_WELL_RE = re.compile(r"^([A-D])([A-H])(0[1-9]|1[0-2])$", re.IGNORECASE)


def plate_letter_from_kit(kit: str) -> str:
    """Return plate letter (A–D) from a TWIST kit name."""
    if kit not in TWIST_KITS:
        raise ValueError(f"Unsupported TWIST kit: {kit!r}")
    return kit.split("-")[1][0].upper()


def normalize_well_id(well: str, kit: str) -> str:
    """
    Normalize a well identifier to standard form (e.g. A01, B03, H12).

    Accepts A1 / A01 / AA01 (Dorado-style, plate letter must match kit).
    """
    raw = well.strip().upper().replace(" ", "")
    if not raw:
        raise ValueError("Empty well ID")

    match = WELL_ID_RE.fullmatch(raw)
    if match:
        row, col = match.group(1).upper(), int(match.group(2))
        return f"{row}{col:02d}"

    dorado = DORADO_WELL_RE.fullmatch(raw)
    if dorado:
        plate, row, col = dorado.group(1).upper(), dorado.group(2).upper(), int(dorado.group(3))
        expected_plate = plate_letter_from_kit(kit)
        if plate != expected_plate:
            raise ValueError(
                f"Well {raw!r} is for plate {plate}, but kit is {kit!r} (plate {expected_plate})"
            )
        return f"{row}{col:02d}"

    raise ValueError(
        f"Unrecognized well ID {well!r}; use e.g. A01 or H12 (or Dorado-style {plate_letter_from_kit(kit)}A01)"
    )


def well_id_to_barcode(well: str, kit: str) -> str:
    """Map a well ID to a Dorado barcode label (barcode01 … barcode96)."""
    normalized = normalize_well_id(well, kit)
    row = normalized[0]
    col = int(normalized[1:])
    row_idx = PLATE_ROWS.index(row)
    index = (col - 1) * 8 + row_idx + 1
    if index < 1 or index > 96:
        raise ValueError(f"Well {well!r} is out of range for a 96-well plate")
    return f"barcode{index:02d}"


def all_well_ids() -> list[str]:
    """Return all 96 well IDs in Dorado column-major order."""
    wells: list[str] = []
    for col in range(1, 13):
        for row in PLATE_ROWS:
            wells.append(f"{row}{col:02d}")
    return wells


def mock_twist_run_rows() -> list[tuple[str, str]]:
    """Mock sample name + well ID pairs for a full 96-well plate."""
    return [(f"Sample{i:02d}", well) for i, well in enumerate(all_well_ids(), start=1)]


def parse_twist_run_file(path: Path, kit: str) -> list[tuple[str, str]]:
    """
    Parse a run file with one sample per line: ``sample,well`` (comma-separated).

    Blank lines and ``#`` comments are ignored.
    """
    rows: list[tuple[str, str]] = []
    seen_samples: set[str] = set()
    seen_wells: set[str] = set()

    for line_no, line in enumerate(path.read_text().splitlines(), start=1):
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        if "," not in stripped:
            raise ValueError(f"{path}:{line_no}: expected 'sample,well' (comma-separated)")
        sample, well = (part.strip() for part in stripped.split(",", 1))
        if not sample:
            raise ValueError(f"{path}:{line_no}: missing sample name")
        if not well:
            raise ValueError(f"{path}:{line_no}: missing well ID")
        if sample in seen_samples:
            raise ValueError(f"{path}:{line_no}: duplicate sample name {sample!r}")
        normalized_well = normalize_well_id(well, kit)
        if normalized_well in seen_wells:
            raise ValueError(f"{path}:{line_no}: duplicate well {normalized_well!r}")
        seen_samples.add(sample)
        seen_wells.add(normalized_well)
        rows.append((sample, normalized_well))

    if not rows:
        raise ValueError(f"No samples found in {path}")
    return rows


def twist_run_to_pod5_rows(kit: str, run_file: Path) -> list[tuple[str, str]]:
    """Convert a TWIST run file to POD5 samplesheet rows (sample, barcode label)."""
    parsed = parse_twist_run_file(run_file, kit)
    return [(sample, well_id_to_barcode(well, kit)) for sample, well in parsed]
