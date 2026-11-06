"""
Build Metatropics samplesheet CSV files (columns: sample, barcode).

Run via the installed console script ``metatropics-samplesheet`` (or ``metatropics_samplesheet``),
or via ``python3 -m metatropics_samplesheet`` (e.g. when running directly from this repo with
``PYTHONPATH=nf-metatropics/assets/submission``).

Modes:
  pod5            — barcode labels for raw POD5 runs (POD5.csv)
  fastq_pass      — barcode labels for basecalled fastq_pass (fastq_pass.csv)
  demultiplexed   — FASTQ paths for demultiplexed reads (samplesheet.csv)
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

from . import __version__
from .twist import TWIST_KITS, mock_twist_run_rows, twist_run_to_pod5_rows

# Longest suffix first so e.g. .fastq.gz is not parsed as .fastq
FASTQ_SUFFIXES = (".fastq.gz", ".fq.gz", ".fastq", ".fq")
BARCODE_PATTERN = re.compile(r"barcode(\d{2})", re.IGNORECASE)

GENERIC_BARCODE_ROWS: list[tuple[str, str]] = [
    ("Sample1", "barcode01"),
    ("Sample2", "barcode02"),
]

HELP_FLAGS = frozenset({"-h", "--help", "-help", "-?"})
KNOWN_MODES = frozenset({"pod5", "fastq_pass", "demultiplexed", "fastq"})
BARCODE_LABEL_MODES = frozenset({"pod5", "fastq_pass"})
DEFAULT_CSV = {
    "pod5": "POD5.csv",
    "fastq_pass": "fastq_pass.csv",
    "demultiplexed": "samplesheet.csv",
    "fastq": "samplesheet.csv",
}


def add_version_arg(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )


def wants_top_level_help(argv: list[str]) -> bool:
    """True for bare help flags; subcommand help (e.g. pod5 -h) stays mode-specific."""
    if not argv:
        return True
    if argv[0] in KNOWN_MODES:
        return False
    return all(arg in HELP_FLAGS for arg in argv)


TOP_LEVEL_EPILOG = """\
examples:
  # TWIST kits: TWIST-96A-UDI, TWIST-96B-UDI, TWIST-96C-UDI, TWIST-96D-UDI
  1. POD5
  # raw POD5 → POD5.csv
  metatropics-samplesheet pod5 -i .
  metatropics-samplesheet pod5 TWIST-96A-UDI
  metatropics-samplesheet pod5 TWIST-96A-UDI wells.csv -o POD5.csv

  2. fastq_pass
  # basecalled fastq_pass → fastq_pass.csv
  metatropics-samplesheet fastq_pass -i .
  metatropics-samplesheet fastq_pass TWIST-96A-UDI
  metatropics-samplesheet fastq_pass TWIST-96A-UDI wells.csv -o fastq_pass.csv

  3. Demultiplexed
  # pass/barcodeNN/ or flat FASTQs → samplesheet.csv
  metatropics-samplesheet demultiplexed -i .
  metatropics-samplesheet demultiplexed TWIST-96A-UDI -o wells.csv
  metatropics-samplesheet demultiplexed -i pass --kit TWIST-96A-UDI --wells wells.csv -o .
"""

FASTQ_EPILOG = TOP_LEVEL_EPILOG

POD5_EPILOG = TOP_LEVEL_EPILOG


def resolve_input_dir(raw: str | None) -> Path:
    path = Path(raw or ".").expanduser()
    if not path.is_dir():
        raise SystemExit(f"Not a directory: {path}")
    return path.resolve()


def resolve_csv_output(
    output: Path | None,
    default_dir: Path,
    default_name: str,
) -> Path:
    """Resolve -o: None or a directory (incl. .) → default_dir/name or DIR/name."""
    if output is None:
        return (default_dir / default_name).resolve()
    path = output.expanduser()
    if str(path) in (".", "..") or path.is_dir():
        return (path / default_name).resolve()
    return path.resolve()


def sanitize_pipeline_sample_id(sample: str) -> str:
    """Strip suffixes the pipeline itself adds to FASTQ/sample names."""
    suffixes = (
        "_fixed",
        "_human_depleted",
        "_host_depleted",
        "_human",
        "_host",
        "_other",
        "_viral",
        "_classification_results",
        ".fastp",
    )
    name = (sample or "").strip()
    while True:
        lowered = name.lower()
        stripped = False
        for suf in suffixes:
            if lowered.endswith(suf.lower()):
                name = name[: -len(suf)]
                stripped = True
                break
        if not stripped:
            break
    name = name.strip("_")
    return name or sample


def sample_name_from_path(path: Path) -> str | None:
    name = path.name
    lower = name.lower()
    for suf in FASTQ_SUFFIXES:
        if lower.endswith(suf):
            return sanitize_pipeline_sample_id(name[: -len(suf)])
    return None


def find_fastq_files(directory: Path) -> list[Path]:
    found: list[Path] = []
    for p in directory.iterdir():
        if not p.is_file():
            continue
        if sample_name_from_path(p) is not None:
            found.append(p.resolve())
    return sorted(found, key=lambda x: x.name.lower())


def find_barcode_folders(directory: Path) -> list[Path]:
    """Return sorted barcodeNN/ subdirectories (MinKNOW pass/ layout)."""
    found: list[Path] = []
    for p in directory.iterdir():
        if p.is_dir() and BARCODE_PATTERN.fullmatch(p.name):
            found.append(p)
    return sorted(found, key=lambda x: int(BARCODE_PATTERN.fullmatch(x.name).group(1)))


def resolve_barcode_folder_input(barcode_dir: Path) -> Path | None:
    """Point at barcodeNN/ if it contains any FASTQ (1 or many; pipeline reads the folder)."""
    files = find_fastq_files(barcode_dir)
    if not files:
        return None
    return barcode_dir.resolve()


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


def write_twist_run_file(
    out_path: Path,
    rows: list[tuple[str, str]],
    *,
    next_hint: str = "the samplesheet CSV",
) -> None:
    """Write sample,well CSV (wells.csv) for editing or direct EPI2ME use with kit_name."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["sample", "well"])
        writer.writerows(rows)
    print(
        f"Wrote {len(rows)} well(s) to {out_path} — edit sample names, "
        f"delete unused rows; use as EPI2ME sample sheet with kit_name, "
        f"or convert to {next_hint}.",
        file=sys.stderr,
    )


def is_twist_kit(name: str) -> bool:
    return name in TWIST_KITS


def run_twist_barcode_csv(
    kit: str,
    run_file: Path | None,
    output: Path | None,
    *,
    default_csv: str,
) -> None:
    if run_file is not None:
        if not run_file.is_file():
            raise SystemExit(f"Not a file: {run_file}")
        try:
            rows = twist_run_to_pod5_rows(kit, run_file)
        except ValueError as exc:
            raise SystemExit(str(exc)) from exc
        out_path = resolve_csv_output(output, Path("."), default_csv)
        write_samplesheet(out_path, rows)
        return

    rows = mock_twist_run_rows()
    out_path = (output or Path("wells.csv")).resolve()
    write_twist_run_file(out_path, rows, next_hint=default_csv)


def build_barcode_path_map(fastq_dir: Path) -> dict[str, Path]:
    """Map barcodeNN → barcode folder path (pipeline reads one or many FASTQs inside)."""
    paths: dict[str, Path] = {}
    for bdir in find_barcode_folders(fastq_dir):
        target = resolve_barcode_folder_input(bdir)
        if target is None:
            print(f"Skipping {bdir.name}: no FASTQ found", file=sys.stderr)
            continue
        paths[bdir.name.lower()] = target
    return paths


def run_fastq_samplesheet(
    fastq_dir: Path,
    output: Path | None,
    *,
    kit: str | None = None,
    wells: Path | None = None,
) -> None:
    """
    Write sample,barcode CSV for demultiplexed reads.

    Prefers MinKNOW/Dorado barcodeNN/ folders (CSV points at the folder).
    Falls back to a flat directory of per-sample FASTQ files.

    With --kit + --wells: lab sample,well sheet → TWIST barcodeNN → folders under pass/.
    """
    out_path = resolve_csv_output(output, fastq_dir, "samplesheet.csv")
    rows: list[tuple[str, str]] = []

    if (kit is None) ^ (wells is None):
        raise SystemExit("Use --kit and --wells together (TWIST plate well → barcode mapping).")

    if kit is not None and wells is not None:
        if not is_twist_kit(kit):
            kits = ", ".join(sorted(TWIST_KITS))
            raise SystemExit(f"Unsupported kit {kit!r}. Expected one of: {kits}")
        if not wells.is_file():
            raise SystemExit(f"Wells file not found: {wells}")
        barcode_paths = build_barcode_path_map(fastq_dir)
        if not barcode_paths:
            raise SystemExit(
                f"No barcodeNN/ folders with FASTQs under {fastq_dir} "
                "(needed to link TWIST wells to demultiplexed reads)"
            )
        try:
            mapped = twist_run_to_pod5_rows(kit, wells)
        except ValueError as exc:
            raise SystemExit(str(exc)) from exc
        skipped: list[tuple[str, str]] = []
        for sample, barcode in mapped:
            path = barcode_paths.get(barcode.lower())
            if path is None:
                skipped.append((sample, barcode))
                continue
            rows.append((sample, str(path)))
        if not rows:
            available = ", ".join(sorted(barcode_paths))
            raise SystemExit(
                f"No wells from {wells.name} match barcode folders under {fastq_dir}. "
                f"Available: {available}"
            )
        if skipped:
            detail = ", ".join(f"{s}→{b}" for s, b in skipped)
            print(
                f"Skipping {len(skipped)} well(s) with no folder under {fastq_dir}: {detail}",
                file=sys.stderr,
            )
        write_samplesheet(out_path, rows)
        print(
            f"Linked {len(rows)} TWIST well(s) ({kit}) to demultiplexed FASTQs under {fastq_dir}",
            file=sys.stderr,
        )
        return

    barcode_paths = build_barcode_path_map(fastq_dir)
    if barcode_paths:
        for name in sorted(barcode_paths, key=lambda x: int(x.replace("barcode", ""))):
            rows.append((name, str(barcode_paths[name])))
        write_samplesheet(out_path, rows)
        return

    paths = find_fastq_files(fastq_dir)
    if not paths:
        raise SystemExit(
            f"No FASTQ files or barcodeNN/ folders found in {fastq_dir} "
            f"(expected *{', *'.join(FASTQ_SUFFIXES)} or barcode01/, …)"
        )

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

    write_samplesheet(out_path, rows)


def run_barcode_label_samplesheet(
    input_dir: Path,
    output: Path | None,
    *,
    default_csv: str,
) -> None:
    """
    Write barcode-label CSV for POD5 or fastq_pass (not file paths).

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
        rows = list(GENERIC_BARCODE_ROWS)
        print(
            f"No barcode labels found in filenames under {input_dir}; "
            f"wrote a generic {default_csv} template — edit sample and barcode columns.",
            file=sys.stderr,
        )

    out_path = resolve_csv_output(output, input_dir, default_csv)
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


def build_demultiplexed_parser(
    subparsers: argparse._SubParsersAction,
) -> argparse.ArgumentParser:
    kits = ", ".join(sorted(TWIST_KITS))
    parser = subparsers.add_parser(
        "demultiplexed",
        aliases=["fastq"],
        help=argparse.SUPPRESS,
        description=(
            "Build samplesheet.csv from demultiplexed FASTQ files.\n"
            "Columns: sample, barcode (full FASTQ path).\n"
            "Optional: --kit + --wells maps TWIST plate well IDs to barcode folders."
        ),
        epilog=FASTQ_EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    add_io_args(
        parser,
        default_help="Directory of demultiplexed FASTQ files or pass/ (default: current directory).",
        output_help="Output CSV or directory (default: DIR/samplesheet.csv; -o . → ./samplesheet.csv).",
    )
    parser.add_argument(
        "--kit",
        default=None,
        metavar="KIT",
        help=(
            f"TWIST UDI kit for well→barcode mapping with --wells "
            f"(one of: {kits})."
        ),
    )
    parser.add_argument(
        "--wells",
        default=None,
        metavar="FILE",
        help="Lab plate file with sample,well lines (e.g. PatientA,A01). Requires --kit.",
    )
    add_version_arg(parser)
    return parser


def build_twist_barcode_parser(mode: str) -> argparse.ArgumentParser:
    kits = ", ".join(sorted(TWIST_KITS))
    default_csv = DEFAULT_CSV[mode]
    parser = argparse.ArgumentParser(
        prog=f"metatropics-samplesheet {mode}",
        description=(
            f"Build a TWIST UDI plate run file or {default_csv} from sample + well IDs.\n"
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
        metavar="wells.csv",
        help="Optional run file (sample,well per line). Without it, a mock template is written.",
    )
    parser.add_argument(
        "-o",
        "--output",
        default=None,
        metavar="FILE",
        help=f"Output file (default: wells.csv or {default_csv}).",
    )
    add_version_arg(parser)
    return parser


def build_barcode_label_parser(
    subparsers: argparse._SubParsersAction,
    mode: str,
) -> argparse.ArgumentParser:
    default_csv = DEFAULT_CSV[mode]
    parser = subparsers.add_parser(
        mode,
        help=argparse.SUPPRESS,
        description=(
            f"Build {default_csv} for {'raw POD5' if mode == 'pod5' else 'basecalled fastq_pass'} data.\n"
            "Columns: sample, barcode (barcode label, e.g. barcode01 — not a file path).\n"
            f"For TWIST UDI plates: {mode} TWIST-96A-UDI [wells.csv]."
        ),
        epilog=POD5_EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    add_io_args(
        parser,
        default_help=f"{'POD5' if mode == 'pod5' else 'fastq_pass'} directory (default: current directory).",
        output_help=f"Output CSV or directory (default: DIR/{default_csv}).",
    )
    add_version_arg(parser)
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
    build_demultiplexed_parser(subparsers)
    build_barcode_label_parser(subparsers, "pod5")
    build_barcode_label_parser(subparsers, "fastq_pass")
    add_version_arg(parser)
    return parser


def write_wells_template(output: Path | None = None) -> None:
    out = (output or Path("wells.csv")).resolve()
    write_twist_run_file(
        out,
        mock_twist_run_rows(),
        next_hint="samplesheet.csv (demultiplexed -i pass/ --kit … --wells wells.csv -o .)",
    )


def dispatch(mode: str, args: argparse.Namespace) -> None:
    input_dir, output = parse_io_args(args)
    if mode in ("demultiplexed", "fastq"):
        wells = Path(args.wells).expanduser() if getattr(args, "wells", None) else None
        run_fastq_samplesheet(
            input_dir,
            output,
            kit=getattr(args, "kit", None),
            wells=wells,
        )
    elif mode in BARCODE_LABEL_MODES:
        run_barcode_label_samplesheet(
            input_dir,
            output,
            default_csv=DEFAULT_CSV[mode],
        )
    else:
        raise SystemExit(f"Unknown mode: {mode}")


def main(argv: list[str] | None = None) -> None:
    argv = list(sys.argv[1:] if argv is None else argv)

    if argv == ["--version"]:
        print(f"metatropics-samplesheet {__version__}")
        return

    parser = build_parser()

    # Top-level help lists modes and examples (-help is treated like --help).
    if wants_top_level_help(argv):
        parser.print_help()
        return

    # TWIST wells.csv: `… demultiplexed TWIST-96A-UDI` or bare `… TWIST-96A-UDI`
    if argv and is_twist_kit(argv[0]):
        out = Path("wells.csv")
        i = 1
        while i < len(argv):
            if argv[i] in HELP_FLAGS:
                parser.print_help()
                return
            if argv[i] in ("-o", "--output") and i + 1 < len(argv):
                out = Path(argv[i + 1]).expanduser()
                i += 2
                continue
            raise SystemExit(
                f"Unexpected argument {argv[i]!r}. "
                "Usage: metatropics-samplesheet TWIST-96A-UDI [-o wells.csv]"
            )
        write_wells_template(out)
        return

    if (
        len(argv) >= 2
        and argv[0] in ("demultiplexed", "fastq")
        and is_twist_kit(argv[1])
    ):
        out = Path("wells.csv")
        rest = argv[2:]
        i = 0
        while i < len(rest):
            if rest[i] in HELP_FLAGS:
                parser.print_help()
                return
            if rest[i] in ("-o", "--output") and i + 1 < len(rest):
                out = Path(rest[i + 1]).expanduser()
                i += 2
                continue
            raise SystemExit(
                f"Unexpected argument {rest[i]!r}. "
                "Usage: metatropics-samplesheet demultiplexed TWIST-96A-UDI [-o wells.csv]"
            )
        write_wells_template(out)
        return

    # TWIST barcode CSV: `… pod5|fastq_pass TWIST-96A-UDI [wells.csv]`
    if len(argv) >= 2 and argv[0] in BARCODE_LABEL_MODES and is_twist_kit(argv[1]):
        mode = argv[0]
        twist_argv = argv[1:]
        twist_parser = build_twist_barcode_parser(mode)
        if twist_argv and twist_argv[0] in HELP_FLAGS:
            twist_parser.print_help()
            return
        args = twist_parser.parse_args(twist_argv)
        if args.kit not in TWIST_KITS:
            twist_parser.error(f"Unsupported kit {args.kit!r}")
        run_file = Path(args.run_file).expanduser() if args.run_file else None
        output = Path(args.output).expanduser() if args.output else None
        run_twist_barcode_csv(
            args.kit,
            run_file,
            output,
            default_csv=DEFAULT_CSV[mode],
        )
        return

    # Bare flags default to demultiplexed: `metatropics-samplesheet -i .`
    if argv[0] not in KNOWN_MODES:
        argv = ["demultiplexed", *argv]

    args = parser.parse_args(argv)
    dispatch(args.command, args)


if __name__ == "__main__":
    main()
