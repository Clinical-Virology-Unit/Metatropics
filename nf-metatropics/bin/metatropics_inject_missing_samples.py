#!/usr/bin/env python3
"""Add samplesheet samples that never reached Virasign to the Metatropics summary.

Empty no-hit rows in CSV + matching HTML dropdown/panels/Sample(s) count.
Independent of which Virasign database was used.
"""
from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from pathlib import Path


def _add_id(ids: list[str], seen: set[str], s: str) -> None:
    s = (s or "").strip()
    if s and s not in seen:
        seen.add(s)
        ids.append(s)


def load_expected_sample_ids(outdir: Path, samplesheet: Path | None) -> list[str]:
    ids: list[str] = []
    seen: set[str] = set()
    candidates: list[Path] = []
    if samplesheet is not None:
        candidates.append(samplesheet)
    candidates.append(outdir / "Summary" / "pipeline_info" / "samplesheet.valid.csv")
    candidates.append(outdir / "samplesheet.valid.csv")
    for sheet in candidates:
        if not sheet.is_file():
            continue
        try:
            with sheet.open("r", newline="", encoding="utf-8") as handle:
                reader = csv.DictReader(handle)
                for row in reader:
                    _add_id(
                        ids,
                        seen,
                        row.get("sample") or row.get("Sample") or row.get("sample_id") or "",
                    )
        except Exception:
            continue
        if ids:
            return ids
    raw_dir = outdir / "Reads" / "raw"
    if raw_dir.is_dir():
        for p in sorted(raw_dir.iterdir()):
            if not p.is_file():
                continue
            name = p.name
            for suf in ("_fixed.fastq.gz", "_fixed.fq.gz", ".fastq.gz", ".fq.gz"):
                if name.endswith(suf):
                    _add_id(ids, seen, name[: -len(suf)])
                    break
    return ids


def _esc(s: str) -> str:
    return (
        str(s)
        .replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


def _replace_json_assign(src: str, const_name: str, new_obj) -> str:
    m = re.search(rf"(?:const|var|let)\s+{re.escape(const_name)}\s*=\s*", src)
    if not m:
        return src
    start = m.end()
    _, end = json.JSONDecoder().raw_decode(src, start)
    return src[:start] + json.dumps(new_obj, separators=(",", ":")) + src[end:]


def inject_csv(summary_csv: Path, expected_ids: list[str]) -> list[str]:
    if not summary_csv.is_file():
        return []
    with summary_csv.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        fields = list(reader.fieldnames or [])
        data = list(reader)
    if not fields:
        return []

    low = {f.lower().strip(): f for f in fields}
    sample_col = next((low[k] for k in ("sample", "sample_id", "sampleid") if k in low), None)
    accession_col = next(
        (low[k] for k in ("accession", "virus", "reference", "reference_id", "ref") if k in low),
        None,
    )
    if not sample_col:
        return []

    present = {(row.get(sample_col) or "").strip() for row in data if (row.get(sample_col) or "").strip()}
    qc_header = "QC reads" if "QC reads" in fields else None
    bg_header = "Background" if "Background" in fields else None

    for sid in expected_ids:
        if sid in present:
            continue
        empty = {f: "" for f in fields}
        empty[sample_col] = sid
        if accession_col:
            empty[accession_col] = ""
        if qc_header:
            empty[qc_header] = "0"
        if bg_header:
            empty[bg_header] = "no"
        data.append(empty)
        present.add(sid)

    with summary_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(data)

    return sorted(present)


def inject_html(
    summary_html: Path,
    expected_ids: list[str],
    db_label: str,
) -> None:
    if not summary_html.is_file() or not expected_ids:
        return
    text = summary_html.read_text(encoding="utf-8")

    m_sd = re.search(r"(?:const|var|let)\s+samplesData\s*=\s*", text)
    m_as = re.search(r"(?:const|var|let)\s+allSamples\s*=\s*", text)
    if not (m_sd and m_as):
        return

    samples_obj, _ = json.JSONDecoder().raw_decode(text, m_sd.end())
    all_obj, _ = json.JSONDecoder().raw_decode(text, m_as.end())
    if not isinstance(samples_obj, dict):
        samples_obj = {}
    if not isinstance(all_obj, list):
        all_obj = []

    meta_template = {
        "database_source": db_label,
        "database_used": "",
        "filtering_criteria": {},
    }
    for dbs in samples_obj.values():
        if not isinstance(dbs, dict):
            continue
        payload = dbs.get(db_label) or next(iter(dbs.values()), None)
        if isinstance(payload, dict) and isinstance(payload.get("metadata"), dict):
            meta_template = dict(payload["metadata"])
            break

    existing_sections = set(re.findall(r'id="sample-([^"]+)"', text))
    existing_options = set(re.findall(r'<option value="([^"]+)">', text))
    to_add_dom = [sid for sid in expected_ids if sid not in existing_sections]

    all_set = set(all_obj)
    for sid in expected_ids:
        if sid not in all_set:
            all_obj.append(sid)
            all_set.add(sid)
        if sid not in samples_obj:
            samples_obj[sid] = {
                db_label: {
                    "references": [],
                    "metadata": meta_template,
                }
            }
    try:
        all_obj = sorted(all_obj, key=lambda x: str(x))
    except Exception:
        pass

    text = _replace_json_assign(text, "samplesData", samples_obj)
    text = _replace_json_assign(text, "allSamples", all_obj)

    n_samples = len(all_obj)
    text, _ = re.subn(
        r'(<div class="number">)\d+(</div>\s*<div class="label">Sample\(s\)</div>)',
        rf"\g<1>{n_samples}\g<2>",
        text,
        count=1,
    )

    if to_add_dom:
        opt_html = "".join(
            f'<option value="{_esc(sid)}">{_esc(sid)} (no viruses detected)</option>\n'
            for sid in to_add_dom
            if sid not in existing_options
        )
        if opt_html:
            m_sel = re.search(
                r'(<select[^>]*id="sampleSelect"[^>]*>)(.*?)(</select>)',
                text,
                flags=re.S,
            )
            if m_sel:
                text = text[: m_sel.start(3)] + opt_html + text[m_sel.start(3) :]

        sections = []
        for sid in to_add_dom:
            es = _esc(sid)
            sections.append(
                f"""
            <div class="sample-section" id="sample-{es}">
                <div class="metadata-section">
                    <h3>📋 Pipeline Configuration</h3>
                    <div class="metadata-grid" id="metadata-{es}">
                        <!-- Metadata will be populated by JavaScript -->
                    </div>
                </div>

                <div style="background: #fff3cd; border: 2px solid #ffc107; border-radius: 8px; padding: 30px; margin: 20px 0; text-align: center;">
                    <h2 style="color: #856404;">⚠️ No Viruses Detected</h2>
                    <p style="color: #856404;">Sample: {es}</p>
                </div>
            </div>
"""
            )
        block = "".join(sections)
        marker = "<!-- Heatmap Section"
        if marker in text:
            text = text.replace(marker, block + "\n            " + marker, 1)
        else:
            text = text.replace(
                '<div class="heatmap-section">',
                block + '\n            <div class="heatmap-section">',
                1,
            )

    summary_html.write_text(text, encoding="utf-8")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--outdir", required=True, type=Path)
    ap.add_argument("--csv", required=True, type=Path)
    ap.add_argument("--html", required=True, type=Path)
    ap.add_argument("--samplesheet", default="", type=str)
    ap.add_argument("--db-label", required=True, type=str)
    args = ap.parse_args()

    sheet = Path(args.samplesheet) if str(args.samplesheet or "").strip() else None
    expected = load_expected_sample_ids(args.outdir, sheet)
    if not expected:
        print("metatropics_inject_missing_samples: no expected sample IDs found; skip", file=sys.stderr)
        return 0

    present = inject_csv(args.csv, expected)
    inject_html(args.html, present or expected, args.db_label)
    print(
        f"metatropics_inject_missing_samples: expected={len(expected)} "
        f"in_report={len(present or expected)} db={args.db_label}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
