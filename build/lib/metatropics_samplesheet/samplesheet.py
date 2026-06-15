"""
Build Metatropics samplesheet CSV files (columns: sample, barcode).

Run via the installed console script ``metatropics-samplesheet`` (or ``metatropics_samplesheet``),
or via ``python3 -m metatropics_samplesheet`` (e.g. when running directly from this repo with
``PYTHONPATH=nf-metatropics/assets/submission``).

Modes:
  fastq (default) — per-sample FASTQ paths (demultiplexed reads)
  pod5            — barcode labels for POD5 or fastq_pass runs (template POD5.csv)
  pod5 TWIST-*    — TWIST UDI plate mapping (well ID → barcode01…96)
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

from .twist import TWIST_KITS, mock_twist_run_rows, twist_run_to_pod5_rows

# Longest suffix first so e.g. .fastq.gz is not parsed as .fastq
FASTQ_SUFFIXES = (".fastq.gz", ".fq.gz", ".fastq", ".fq")
BARCODE_PATTERN = re.compile(r"barcode(\d{2})", re.IGNORECASE)

GENERIC_POD5_ROWS: list[tuple[str, str]] = [
    ("Sample1", "barcode01"),
    ("Sample2", "barcode02"),
]

HELP_FLAGS = frozenset({"-h", "--help", "-help", "-?"})


def wants_top_level_help(argv: list[str]) -> bool:
    """True for bare help flags; subcommand help (e.g. pod5 -h) stays mode-specific."""
    if not argv:
        return True
    if argv[0] in ("fastq", "pod5"):
        return False
    return all(arg in HELP_FLAGS for arg in argv)


TOP_LEVEL_EPILOG = """\
examples:
  metatropics-samplesheet -i .
  metatropics-samplesheet pod5 -i .
  metatropics-samplesheet pod5 TWIST-96A-UDI
  metatropics-samplesheet pod5 TWIST-96A-UDI run.txt -o POD5.csv
"""

FASTQ_EPILOG = TOP_LEVEL_EPILOG + """\

fastq mode writes samplesheet.csv with one row per FASTQ file in the directory (top level only).
use the output CSV as --input in params_fastq.yaml.
"""

POD5_EPILOG = TOP_LEVEL_EPILOG + """\

pod5 mode writes POD5.csv with barcode labels (not file paths). edit sample names and barcodes,
then use the output CSV as --input in params_POD5.yaml or params_fastq_pass.yaml.

TWIST UDI plates (TWIST-96A-UDI … TWIST-96D-UDI):
  metatropics-samplesheet pod5 TWIST-96A-UDI
      → writes run.txt (mock Sample01,A01 … Sample96,H12)
  metatropics-samplesheet pod5 TWIST-96A-UDI run.txt
      → writes POD5.csv from comma-separated sample,well lines
  Set the same kit in params (kit_name) when running the pipeline.
"""


def resolve_input_dir(raw: str | None) -> Path:
    path = Path(raw or ".").expanduser()
    if not path.is_dir():
        raise SystemExit(f"Not a directory: {path}")
    return path.resolve()


def sample_name_from_path(path: Path) -> str | None:
    name = path.name
    lower = name.lower()
    for suf in FASTQ_SUFFIXES:
        if lower.endswith(suf):
            return name[: -len(suf)]
    return None


def find_fastq_files(directory: Path) -> list[Path]:
    found: list[Path] = []
    for p in directory.iterdir():
        if not p.is_file():
            continue
        if sample_name_from_path(p) is not None:
            found.append(p.resolve())
    return sorted(found, key=lambda x: x.name.lower())


def find_barcodes_in_dir(directory: Path) -> list[str]:
    """Return sorted barcode labels (e.g. barcode01) found in filenames."""
    found: set[str] = set()
    for p in directory.iterdir():
        if not p.is_file():
            continue
        for match in BARCODE_PATTERN.finditer(p.name):
            found.add(f"barcode{match.group(1)}")
    return sorted(found, key=lambda x: int(x.replace("barcode", "")))


def write_samplesheet(out_path: Path, rows: list[tuple[str, str]]) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["sample", "barcode"])
        writer.writerows(rows)
    print(f"Wrote {len(rows)} sample(s) to {out_path}", file=sys.stderr)


def write_twist_run_file(out_path: Path, rows: list[tuple[str, str]]) -> None:
    """Write sample,well text for manual editing before POD5.csv generation."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["# sample,well  (edit sample names; remove unused wells)", ""]
    lines.extend(f"{sample},{well}" for sample, well in rows)
    out_path.write_text("\n".join(lines) + "\n")
    print(
        f"Wrote {len(rows)} well(s) to {out_path} — edit sample names, "
        "delete unused rows, then run with this file to create POD5.csv.",
        file=sys.stderr,
    )


def is_twist_kit(name: str) -> bool:
    return name in TWIST_KITS


def run_twist_pod5(kit: str, run_file: Path | None, output: Path | None) -> None:
    if run_file is not None:
        if not run_file.is_file():
            raise SystemExit(f"Not a file: {run_file}")
        try:
            rows = twist_run_to_pod5_rows(kit, run_file)
        except ValueError as exc:
            raise SystemExit(str(exc)) from exc
        out_path = (output or Path("POD5.csv")).resolve()
        write_samplesheet(out_path, rows)
        return

    rows = mock_twist_run_rows()
    out_path = (output or Path("run.txt")).resolve()
    write_twist_run_file(out_path, rows)


def run_fastq_samplesheet(fastq_dir: Path, output: Path | None) -> None:
    """Write sample,barcode CSV for one directory of FASTQ files (top level only)."""
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

    out_path = (output or fastq_dir / "samplesheet.csv").resolve()
    write_samplesheet(out_path, rows)


def run_pod5_samplesheet(input_dir: Path, output: Path | None) -> None:
    """
    Write POD5.csv for raw POD5 or fastq_pass runs (barcode labels, not file paths).

    If filenames contain barcodeNN labels, one row is written per detected barcode.
    Otherwise a small generic template is written for manual editing.
    """
    barcodes = find_barcodes_in_dir(input_dir)
    if barcodes:
        rows = [(f"Sample{i}", barcode) for i, barcode in enumerate(barcodes, start=1)]
        print(
            f"Detected {len(barcodes)} barcode label(s) in {input_dir}; "
            "edit sample names in the CSV as needed.",
            file=sys.stderr,
        )
    else:
        rows = list(GENERIC_POD5_ROWS)
        print(
            f"No barcode labels found in filenames under {input_dir}; "
            "wrote a generic POD5.csv template — edit sample and barcode columns.",
            file=sys.stderr,
        )

    out_path = (output or input_dir / "POD5.csv").resolve()
    write_samplesheet(out_path, rows)


def add_io_args(
    parser: argparse.ArgumentParser,
    *,
    default_help: str,
    output_help: str,
) -> None:
    parser.add_argument(
        "-i",
        "--input",
        default=None,
        metavar="DIR",
        help=default_help,
    )
    parser.add_argument(
        "-o",
        "--output",
        default=None,
        metavar="FILE",
        help=output_help,
    )
    parser.add_argument(
        "dir_pos",
        nargs="?",
        default=None,
        metavar="DIR",
        help="Input directory (overrides -i if given).",
    )


def parse_io_args(args: argparse.Namespace) -> tuple[Path, Path | None]:
    raw = (
        args.dir_pos
        if args.dir_pos is not None
        else (args.input if args.input is not None else ".")
    )
    input_dir = resolve_input_dir(raw)
    output = Path(args.output).expanduser() if args.output else None
    return input_dir, output


def build_fastq_parser(subparsers: argparse._SubParsersAction) -> argparse.ArgumentParser:
    parser = subparsers.add_parser(
        "fastq",
        help=argparse.SUPPRESS,
        description=(
            "Build a Metatropics samplesheet from per-sample demultiplexed FASTQ files.\n"
            "Columns: sample, barcode (full FASTQ path)."
        ),
        epilog=FASTQ_EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    add_io_args(
        parser,
        default_help="Directory of demultiplexed FASTQ files (default: current directory).",
        output_help="Output CSV (default: DIR/samplesheet.csv).",
    )
    return parser


def build_twist_pod5_parser() -> argparse.ArgumentParser:
    kits = ", ".join(sorted(TWIST_KITS))
    parser = argparse.ArgumentParser(
        prog="metatropics-samplesheet pod5",
        description=(
            "Build a TWIST UDI plate run file or POD5.csv from sample + well IDs.\n"
            f"Supported kits: {kits}."
        ),
        epilog=POD5_EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "kit",
        metavar="KIT",
        help=f"TWIST kit name (e.g. TWIST-96A-UDI). One of: {kits}.",
    )
    parser.add_argument(
        "run_file",
        nargs="?",
        default=None,
        metavar="RUN.txt",
        help="Optional run file (sample,well per line). Without it, a mock template is written.",
    )
    parser.add_argument(
        "-o",
        "--output",
        default=None,
        metavar="FILE",
        help="Output file (default: run.txt or POD5.csv).",
    )
    return parser


def build_pod5_parser(subparsers: argparse._SubParsersAction) -> argparse.ArgumentParser:
    parser = subparsers.add_parser(
        "pod5",
        help=argparse.SUPPRESS,
        description=(
            "Build a Metatropics POD5.csv template for raw POD5 or basecalled fastq_pass data.\n"
            "Columns: sample, barcode (barcode label, e.g. barcode01 — not a file path).\n"
            "For TWIST UDI plates, pass a kit name: pod5 TWIST-96A-UDI [run.txt]."
        ),
        epilog=POD5_EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    add_io_args(
        parser,
        default_help="POD5 or fastq_pass directory (default: current directory).",
        output_help="Output CSV (default: DIR/POD5.csv).",
    )
    return parser


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="metatropics-samplesheet",
        description=(
            "Build Metatropics samplesheet CSV files (columns: sample, barcode)."
        ),
        epilog=TOP_LEVEL_EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    subparsers = parser.add_subparsers(dest="command", help=argparse.SUPPRESS)
    build_fastq_parser(subparsers)
    build_pod5_parser(subparsers)
    return parser


def dispatch(mode: str, args: argparse.Namespace) -> None:
    input_dir, output = parse_io_args(args)
    if mode == "fastq":
        run_fastq_samplesheet(input_dir, output)
    elif mode == "pod5":
        run_pod5_samplesheet(input_dir, output)
    else:
        raise SystemExit(f"Unknown mode: {mode}")


def main(argv: list[str] | None = None) -> None:
    argv = list(sys.argv[1:] if argv is None else argv)
    parser = build_parser()

    # Top-level help lists both modes and examples (-help is treated like --help).
    if wants_top_level_help(argv):
        parser.print_help()
        return

    # TWIST pod5: `metatropics-samplesheet pod5 TWIST-96A-UDI [run.txt]`
    if len(argv) >= 2 and argv[0] == "pod5" and is_twist_kit(argv[1]):
        twist_argv = argv[1:]
        twist_parser = build_twist_pod5_parser()
        if twist_argv and twist_argv[0] in HELP_FLAGS:
            twist_parser.print_help()
            return
        args = twist_parser.parse_args(twist_argv)
        if args.kit not in TWIST_KITS:
            twist_parser.error(f"Unsupported kit {args.kit!r}")
        run_file = Path(args.run_file).expanduser() if args.run_file else None
        output = Path(args.output).expanduser() if args.output else None
        run_twist_pod5(args.kit, run_file, output)
        return

    # Backward compatible: `metatropics-samplesheet -i .` runs fastq mode.
    if argv[0] not in ("fastq", "pod5"):
        argv = ["fastq", *argv]

    args = parser.parse_args(argv)
    dispatch(args.command, args)


if __name__ == "__main__":
    main()
