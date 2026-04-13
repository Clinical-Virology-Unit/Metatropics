"""
Build a Metatropics FASTQ samplesheet (columns: sample, barcode).

Console entry point: ``metatropics-samplesheet`` (after ``pip install -e .``).
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

# Longest suffix first so e.g. .fastq.gz is not parsed as .fastq
FASTQ_SUFFIXES = (".fastq.gz", ".fq.gz", ".fastq", ".fq")


def sample_name_from_path(path: Path) -> str | None:
    name = path.name
    lower = name.lower()
    for suf in FASTQ_SUFFIXES:
        if lower.endswith(suf):
            return name[: -len(suf)]
    return None


def find_fastq_files(directory: Path) -> list[Path]:
    if not directory.is_dir():
        raise SystemExit(f"Not a directory: {directory}")
    found: list[Path] = []
    for p in directory.iterdir():
        if not p.is_file():
            continue
        if sample_name_from_path(p) is not None:
            found.append(p.resolve())
    return sorted(found, key=lambda x: x.name.lower())


def run_samplesheet(fastq_dir: Path, output: Path | None) -> None:
    """Write sample,barcode CSV for one directory of FASTQ files (top level only)."""
    fastq_dir = fastq_dir.expanduser().resolve()
    paths = find_fastq_files(fastq_dir)
    if not paths:
        raise SystemExit(
            f"No FASTQ files found in {fastq_dir} (expected *{', *'.join(FASTQ_SUFFIXES)})"
        )

    rows: list[tuple[str, str]] = []
    seen_samples: dict[str, Path] = {}
    for fp in paths:
        stem = sample_name_from_path(fp)
        assert stem is not None
        if stem in seen_samples:
            raise SystemExit(
                f"Duplicate sample name {stem!r} from:\n  {seen_samples[stem]}\n  {fp}"
            )
        seen_samples[stem] = fp
        rows.append((stem, str(fp)))

    if output is None:
        out_path = fastq_dir / "samplesheet.csv"
    else:
        out_path = output.expanduser().resolve()
        out_path.parent.mkdir(parents=True, exist_ok=True)

    with out_path.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["sample", "barcode"])
        w.writerows(rows)

    print(f"Wrote {len(rows)} sample(s) to {out_path}", file=sys.stderr)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build Metatropics FASTQ samplesheet (columns: sample, barcode).",
    )
    parser.add_argument(
        "-i",
        "--input",
        default=None,
        metavar="DIR",
        help="Directory of FASTQ files (default: current directory).",
    )
    parser.add_argument(
        "-o",
        "--output",
        default=None,
        metavar="FILE",
        help="Output CSV (default: DIR/samplesheet.csv).",
    )
    parser.add_argument(
        "dir_pos",
        nargs="?",
        default=None,
        metavar="DIR",
        help="FASTQ directory (overrides -i if given).",
    )
    args = parser.parse_args()

    raw = (
        args.dir_pos
        if args.dir_pos is not None
        else (args.input if args.input is not None else ".")
    )
    input_path = Path(raw).expanduser()
    if not input_path.is_dir():
        sys.exit(f"Not a directory: {input_path}")

    out = Path(args.output) if args.output else None
    run_samplesheet(input_path, out)


if __name__ == "__main__":
    main()
