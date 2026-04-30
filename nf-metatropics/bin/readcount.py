#!/usr/bin/env python3
"""
Compute read-category counts (raw/trimmed/human/host/viral/non-viral) and render:
  - read_counts.csv
  - read_distribution.html (interactive: read mix + Virasign viral-by-taxon breakdown)
  - read_distribution.pdf  (static export)

Viral reads are computed from Virasign *confident hits only*:
  Classification/virasign/*/*/*_final_selected_references.json

This script is designed to be run from a working directory that contains
the `read_count/` staging folder.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Tuple

# Stacked category order (bottom → top in horizontal bars): broad background → focal viral
CATEGORY_ORDER = ["trimmed_reads", "human_reads", "host_reads", "non_viral", "viral"]
CATEGORY_LABELS = {
    "trimmed_reads": "Quality trimmed",
    "human_reads": "Human host",
    "host_reads": "Other hosts",
    "non_viral": "Non-viral",
    "viral": "Viral",
}
# Distinct, print-safe palette (OK for colour-blind workflows when not sole cue)
CATEGORY_COLORS = {
    "trimmed_reads": "#94a3b8",
    "human_reads": "#fb923c",
    "host_reads": "#38bdf8",
    "non_viral": "#cbd5e1",
    "viral": "#c026d3",
}


def strip_extensions(name: str) -> str:
    while True:
        new = re.sub(r"\.(fastq|fastp|fq|gz|meta|csv)$", "", name)
        if new == name:
            return name
        name = new


def extract_sample_name(filename: str) -> str:
    name = strip_extensions(filename)
    # New depletion suffixes (published by SAMTOOLS_hFASTQ / SAMTOOLS_hoFASTQ).
    name = re.sub(r"_(human|host)_depleted$", "", name)
    name = re.sub(r"_(human|host)$", "", name)
    # Legacy suffix kept for backwards compatibility.
    name = re.sub(r"_other$", "", name)
    name = re.sub(r"_viral$", "", name)
    name = re.sub(r"_classification_results$", "", name)
    name = re.sub(r"_fixed$", "", name)
    name = re.sub(r"_T\d+$", "", name)
    name = re.sub(r"\.fastp$", "", name)
    return name


def read_samplesheet_order(outdir: Path) -> List[str]:
    """
    Return sample names in the same order as the (validated) submission samplesheet.

    We prefer the validated samplesheet because it is the canonical pipeline input and
    matches what `SAMPLESHEET_CHECK_METATROPICS` produced.
    """
    candidates = [
        outdir / "Summary" / "pipeline_info" / "samplesheet.valid.csv",
        outdir / "samplesheet.valid.csv",
    ]
    sheet = next((p for p in candidates if p.exists()), None)
    if sheet is None:
        return []

    order: List[str] = []
    try:
        with sheet.open(newline="") as fh:
            r = csv.DictReader(fh)
            if not r.fieldnames or "sample" not in r.fieldnames:
                return []
            for row in r:
                raw = (row.get("sample") or "").strip()
                if not raw:
                    continue
                # Normalize to the same "logical sample" key used elsewhere in this script.
                # This keeps ordering stable even if input names include legacy suffixes.
                order.append(extract_sample_name(raw))
    except Exception:
        return []

    # Preserve first occurrence order while de-duplicating.
    seen = set()
    out: List[str] = []
    for s in order:
        if s in seen:
            continue
        seen.add(s)
        out.append(s)
    return out


def count_fastq_gz_reads(path: Path) -> int:
    n = 0
    try:
        with gzip.open(path, "rt", errors="replace") as fh:
            while True:
                l1 = fh.readline()
                if not l1:
                    break
                fh.readline()
                fh.readline()
                fh.readline()
                n += 1
    except EOFError as e:
        # This almost always indicates a truncated/corrupt .gz (incomplete copy/download).
        raise RuntimeError(f"Corrupt gzip FASTQ (unexpected EOF): {path}") from e
    return n


def count_dir(pattern: re.Pattern[str], directory: Path) -> Dict[str, int]:
    out: Dict[str, int] = {}
    if not directory.exists():
        return out
    for p in directory.iterdir():
        if not p.is_file():
            continue
        if not pattern.search(p.name):
            continue
        sample = extract_sample_name(p.name)
        out[sample] = out.get(sample, 0) + count_fastq_gz_reads(p)
    return out


def _virus_label_from_record(rec: Dict[str, Any]) -> str:
    org = rec.get("organism")
    if org and str(org).strip():
        return str(org).strip()
    sp = rec.get("viral_species")
    if sp and str(sp).strip():
        return str(sp).strip()
    acc = rec.get("accession")
    if acc and str(acc).strip():
        return str(acc).strip()
    desc = rec.get("description")
    if desc and str(desc).strip():
        s = str(desc).strip()
        return (s[:72] + "…") if len(s) > 75 else s
    return "Unknown"


def parse_viral_species_by_sample(outdir: Path) -> Dict[str, List[Tuple[str, int]]]:
    """Per sample: list of (virus_label, mapped_reads) from *_final_selected_references.json (confident only)."""
    root = outdir / "Classification" / "virasign"
    if not root.exists():
        return {}
    files = list(root.glob("*/*/*_final_selected_references.json"))
    acc: Dict[str, Dict[str, int]] = {}
    for f in sorted(files):
        raw_name = f.name.replace("_final_selected_references.json", "")
        sample = extract_sample_name(raw_name)
        try:
            data = json.loads(f.read_text())
        except Exception:
            data = None
        if not isinstance(data, list):
            continue
        bucket = acc.setdefault(sample, {})
        for rec in data:
            if not isinstance(rec, dict):
                continue
            try:
                mr = int(rec.get("mapped_reads", 0))
            except Exception:
                mr = 0
            if mr <= 0:
                continue
            lab = _virus_label_from_record(rec)
            bucket[lab] = bucket.get(lab, 0) + mr
    out: Dict[str, List[Tuple[str, int]]] = {}
    for sample, labmap in acc.items():
        out[sample] = sorted(labmap.items(), key=lambda x: -x[1])
    return out


def parse_viral_reads_from_virasign(outdir: Path) -> Dict[str, int]:
    root = outdir / "Classification" / "virasign"
    if not root.exists():
        return {}
    files = list(root.glob("*/*/*_final_selected_references.json"))
    out: Dict[str, int] = {}
    for f in sorted(files):
        sample = f.name.replace("_final_selected_references.json", "")
        total = 0
        try:
            data = json.loads(f.read_text())
        except Exception:
            data = None
        if isinstance(data, list):
            for rec in data:
                if isinstance(rec, dict):
                    v = rec.get("mapped_reads")
                    try:
                        total += int(v)
                    except Exception:
                        pass
        key = extract_sample_name(sample)
        out[key] = out.get(key, 0) + total
    return out


def safe_get(d: Dict[str, int], k: str) -> int:
    return int(d.get(k, 0))


@dataclass(frozen=True)
class Row:
    sample: str
    raw: int
    trimmed: int
    human_depleted: int
    host_depleted: int
    viral: int


def compute_rows(
    read_count_dir: Path,
    outdir: Path,
    host_status: str,
) -> Tuple[List[Row], bool, bool]:
    # Raw reads can be named either `*_fixed.fastq.gz` (older convention) or `*.fastq.gz`.
    # Exclude fastp outputs which are counted separately.
    raw = count_dir(re.compile(r"(?<!\\.fastp)(?:_fixed)?\\.fastq\\.gz$"), read_count_dir)
    trimmed = count_dir(re.compile(r"\.fastp\.fastq\.gz$"), read_count_dir)
    human_dep = count_dir(re.compile(r"\.fastq\.gz$"), read_count_dir / "nohuman")
    host_dep = count_dir(re.compile(r"\.fastq\.gz$"), read_count_dir / "nohost")
    viral = parse_viral_reads_from_virasign(outdir)

    # If the staging folder is empty (common when manually testing), fall back to counting
    # directly from the pipeline outdir's Reads/* structure.
    if not raw and not trimmed:
        reads_root = outdir / "Reads"
        raw = count_dir(re.compile(r"(?<!\\.fastp)(?:_fixed)?\\.fastq\\.gz$"), reads_root / "fix")
        trimmed = count_dir(re.compile(r"\.fastp\.fastq\.gz$"), reads_root / "fastplong")
        # Depletion folders now publish *_human_depleted.fastq.gz / *_host_depleted.fastq.gz.
        # Keep the legacy *_other.fastq.gz pattern so older outdirs still count correctly.
        human_dep = count_dir(re.compile(r"_(human_depleted|other)\.fastq\.gz$"), reads_root / "nohuman")
        host_dep = count_dir(re.compile(r"_(host_depleted|other)\.fastq\.gz$"), reads_root / "nohost")

    present = {*raw.keys(), *trimmed.keys(), *human_dep.keys(), *host_dep.keys(), *viral.keys()}
    order = read_samplesheet_order(outdir)
    samples = [s for s in order if s in present]
    # Append any samples not in the samplesheet (e.g. manual tests / extra files) deterministically.
    extras = sorted(present.difference(samples))
    samples.extend(extras)

    include_human = host_status in {"human_only", "both"} and len(human_dep) > 0
    include_other = host_status in {"other_only", "both"} and len(host_dep) > 0

    rows: List[Row] = []
    for s in samples:
        r = safe_get(raw, s)
        t = safe_get(trimmed, s)
        hd = safe_get(human_dep, s)
        od = safe_get(host_dep, s)
        v = safe_get(viral, s)

        if not include_human:
            hd = t
        if not include_other:
            od = hd

        rows.append(Row(sample=s, raw=r, trimmed=t, human_depleted=hd, host_depleted=od, viral=v))
    return rows, include_human, include_other


def write_csv(rows: List[Row], include_human: bool, include_other: bool, out_csv: Path) -> None:
    header = [
        "sample",
        "raw",
        "trimmed_reads",
        "trimmed_reads_pct",
    ]
    if include_human:
        header += ["human_reads", "human_reads_pct"]
    if include_other:
        header += ["host_reads", "host_reads_pct"]
    header += ["viral", "viral_pct", "non_viral", "non_viral_pct"]

    def denom_for(row: Row) -> int:
        # Prefer true raw reads, but keep the plot usable in partial test outputs.
        if row.raw > 0:
            return row.raw
        return max(row.trimmed, row.human_depleted, row.host_depleted, row.viral, 1)

    def pct(num: int, denom: int) -> float:
        return round((num / denom * 100.0) if denom > 0 else 0.0, 2)

    with out_csv.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        for row in rows:
            denom = denom_for(row)
            trimmed_reads = max(row.raw - row.trimmed, 0)
            human_reads = max(row.trimmed - row.human_depleted, 0) if include_human else 0
            host_reads = max(row.human_depleted - row.host_depleted, 0) if include_other else 0
            non_viral = max(
                (row.host_depleted if include_other else row.human_depleted if include_human else row.trimmed) - row.viral,
                0,
            )

            rec = [
                row.sample,
                row.raw,
                trimmed_reads,
                pct(trimmed_reads, denom),
            ]
            if include_human:
                rec += [human_reads, pct(human_reads, denom)]
            if include_other:
                rec += [host_reads, pct(host_reads, denom)]
            rec += [
                row.viral,
                pct(row.viral, denom),
                non_viral,
                pct(non_viral, denom),
            ]
            w.writerow(rec)


def _stacked_read_figure(
    df: Any,
    x_cols: List[str],
    *,
    x_is_percent: bool,
    height_px: int,
) -> Any:
    import plotly.graph_objects as go  # type: ignore

    samples = df["sample"].astype(str).tolist()
    traces: List[Any] = []
    for col in x_cols:
        key = col.replace("_pct", "") if str(col).endswith("_pct") else str(col)
        color = CATEGORY_COLORS.get(key, "#64748b")
        name = CATEGORY_LABELS.get(key, key)
        if x_is_percent:
            xvals = df[col].fillna(0).astype(float).tolist()
            hover_tmpl = "<b>%{y}</b><br>" + name + ": %{x:.1f}%<extra></extra>"
        else:
            xvals = [int(round(float(v))) for v in df[col].fillna(0).tolist()]
            hover_tmpl = "<b>%{y}</b><br>" + name + ": %{x:.0f} reads<extra></extra>"
        traces.append(
            go.Bar(
                name=name,
                x=xvals,
                y=samples,
                orientation="h",
                marker=dict(color=color, line=dict(width=0)),
                hovertemplate=hover_tmpl,
            )
        )

    x_title = "Percentage (%)" if x_is_percent else "Read count"
    fig = go.Figure(data=traces)
    # Title is rendered in HTML above the chart (light theme); keep plot area uncluttered.
    fig.update_layout(
        title=dict(text=""),
        annotations=[],
        barmode="stack",
        bargap=0.35,
        height=height_px,
        dragmode=False,
        paper_bgcolor="#ffffff",
        plot_bgcolor="#f8fafc",
        font=dict(color="#0f172a", family="system-ui, Segoe UI, Roboto, Helvetica, Arial, sans-serif", size=12),
        xaxis=dict(
            title=dict(text=x_title, font=dict(color="#475569", size=12)),
            tickfont=dict(color="#334155"),
            gridcolor="#e2e8f0",
            zeroline=False,
            ticksuffix="" if x_is_percent else "",
            tickformat=".0f" if x_is_percent else ",",
            separatethousands=True,
            range=[0, 100] if x_is_percent else None,
        ),
        yaxis=dict(
            title="",
            automargin=True,
            categoryorder="array",
            categoryarray=samples,
            tickfont=dict(color="#334155"),
        ),
        legend=dict(
            font=dict(color="#334155", size=11),
            bgcolor="rgba(255,255,255,0.96)",
            bordercolor="#e2e8f0",
            borderwidth=1,
            orientation="h",
            yanchor="top",
            y=-0.17,
            xanchor="center",
            x=0.5,
            traceorder="normal",
        ),
        margin=dict(l=120, r=28, t=36, b=132),
        hovermode="closest",
    )
    return fig


def _viral_detail_figure(
    per_sample: Dict[str, List[Tuple[str, int]]],
    sample_order: List[str],
    height_px: int,
    *,
    as_fraction: bool,
) -> Any:
    """Stacked horizontal bars by taxon: absolute reads, or 100% within-sample (Virasign confident JSON)."""
    import pandas as pd  # type: ignore
    import plotly.express as px  # type: ignore
    import plotly.graph_objects as go  # type: ignore

    max_taxa = 22
    rows: List[Dict[str, Any]] = []
    for sample in sample_order:
        pairs = list(per_sample.get(sample, []))
        if not pairs:
            continue
        pairs = sorted(pairs, key=lambda x: -x[1])
        if len(pairs) > max_taxa:
            head = pairs[: max_taxa - 1]
            other_reads = sum(r for _, r in pairs[max_taxa - 1 :])
            if other_reads > 0:
                pairs = head + [("Other viruses", other_reads)]
            else:
                pairs = head
        for lab, reads in pairs:
            if reads > 0:
                rows.append({"sample": sample, "virus": lab, "reads": reads})

    if not rows:
        fig = go.Figure()
        fig.add_annotation(
            text="No confident viral hits (*_final_selected_references.json)",
            xref="paper",
            yref="paper",
            x=0.5,
            y=0.55,
            showarrow=False,
            font=dict(size=14, color="#64748b"),
        )
        fig.update_layout(
            paper_bgcolor="#ffffff",
            plot_bgcolor="#f8fafc",
            font=dict(color="#0f172a", family="system-ui, Segoe UI, Roboto, Helvetica, Arial, sans-serif", size=12),
            height=max(380, min(height_px, 520)),
            margin=dict(l=80, r=80, t=36, b=60),
            xaxis=dict(visible=False),
            yaxis=dict(visible=False),
            dragmode=False,
        )
        return fig

    df = pd.DataFrame(rows)
    vorder = df.groupby("virus")["reads"].sum().sort_values(ascending=False).index.tolist()
    tot = df.groupby("sample")["reads"].transform("sum")
    df["pct_viral"] = (df["reads"] / tot.replace(0, float("nan")) * 100.0).fillna(0.0).round(4)

    xcol = "pct_viral" if as_fraction else "reads"
    fig = px.bar(
        df,
        x=xcol,
        y="sample",
        color="virus",
        orientation="h",
        category_orders={"sample": list(sample_order), "virus": list(vorder)},
        color_discrete_sequence=(px.colors.qualitative.Bold + px.colors.qualitative.Dark24) * 4,
        labels={"sample": "Sample", "virus": ""},
    )
    for tr in fig.data:
        vname = str(getattr(tr, "name", "") or "")
        ys = list(tr.y) if getattr(tr, "y", None) is not None else []
        cds: List[Any] = []
        for yb in ys:
            m = df[(df["virus"] == vname) & (df["sample"] == yb)]
            if len(m):
                r = int(m["reads"].iloc[0])
                p = float(m["pct_viral"].iloc[0])
            else:
                r, p = 0, 0.0
            cds.append([r] if as_fraction else [p])
        if as_fraction:
            tr.update(
                hovertemplate="<b>%{y}</b><br>"
                + vname
                + ": %{x:.2f}% of viral<br>%{customdata[0]:,} reads<extra></extra>",
                customdata=cds,
            )
        else:
            tr.update(
                hovertemplate="<b>%{y}</b><br>"
                + vname
                + ": %{x:,} reads<br>%{customdata[0]:.1f}% of viral in sample<extra></extra>",
                customdata=cds,
            )
    nh = int(max(420, min(2400, len(sample_order) * 26 + 260)))
    xaxis_reads = dict(
        title="Viral reads",
        title_font=dict(size=12, color="#475569"),
        tickfont=dict(color="#334155"),
        gridcolor="#e2e8f0",
        separatethousands=True,
    )
    xaxis_frac = dict(
        title="% of viral reads (per sample)",
        title_font=dict(size=12, color="#475569"),
        tickfont=dict(color="#334155"),
        gridcolor="#e2e8f0",
        range=[0, 100],
        ticksuffix="",
        tickformat=".0f",
    )
    fig.update_layout(
        barmode="stack",
        bargap=0.32,
        paper_bgcolor="#ffffff",
        plot_bgcolor="#f8fafc",
        font=dict(color="#0f172a", family="system-ui, Segoe UI, Roboto, Helvetica, Arial, sans-serif", size=11),
        xaxis=xaxis_frac if as_fraction else xaxis_reads,
        yaxis=dict(title="", tickfont=dict(color="#334155")),
        legend=dict(
            title=dict(text=""),
            orientation="v",
            yanchor="middle",
            y=0.5,
            x=1.01,
            xanchor="left",
            bgcolor="rgba(255,255,255,0.96)",
            bordercolor="#e2e8f0",
            borderwidth=1,
            font=dict(size=9),
            traceorder="normal",
        ),
        margin=dict(l=120, r=220, t=36, b=120),
        height=nh,
        dragmode=False,
        hovermode="closest",
    )
    return fig


def render_plotly_html(html_out: Path, csv_path: Path, outdir: Path) -> None:
    try:
        import pandas as pd  # type: ignore
    except Exception as e:
        sys.stderr.write(f"WARNING: pandas not available ({e}). Skipping HTML rendering.\n")
        return
    try:
        import plotly.graph_objects as go  # type: ignore  # noqa: F401
    except Exception as e:
        sys.stderr.write(f"WARNING: Plotly not available ({e}). Skipping HTML rendering.\n")
        return

    df = pd.read_csv(csv_path)
    pct_cols = [f"{c}_pct" for c in CATEGORY_ORDER if f"{c}_pct" in df.columns]
    count_cols = [c for c in CATEGORY_ORDER if c in df.columns]
    if not pct_cols:
        sys.stderr.write("WARNING: No percentage columns in CSV; skipping HTML.\n")
        return

    n = max(len(df), 1)
    height_px = int(max(460, min(2400, 28 * n + 220)))

    fig_pct = _stacked_read_figure(df, pct_cols, x_is_percent=True, height_px=height_px)
    fig_cnt = _stacked_read_figure(df, count_cols, x_is_percent=False, height_px=height_px)

    from plotly.utils import PlotlyJSONEncoder  # type: ignore

    spec_pct = fig_pct.to_plotly_json()
    spec_cnt = fig_cnt.to_plotly_json()
    for _spec in (spec_pct, spec_cnt):
        _lay = _spec.setdefault("layout", {})
        _lay.pop("config", None)
        _lay["dragmode"] = False
    pct_json = json.dumps(spec_pct, cls=PlotlyJSONEncoder, allow_nan=False)
    cnt_json = json.dumps(spec_cnt, cls=PlotlyJSONEncoder, allow_nan=False)

    viral_by_sample = parse_viral_species_by_sample(outdir)
    sample_order = df["sample"].astype(str).tolist()
    fig_vir_reads = _viral_detail_figure(viral_by_sample, sample_order, height_px, as_fraction=False)
    fig_vir_frac = _viral_detail_figure(viral_by_sample, sample_order, height_px, as_fraction=True)
    spec_vir_reads = fig_vir_reads.to_plotly_json()
    spec_vir_frac = fig_vir_frac.to_plotly_json()
    for _spec in (spec_vir_reads, spec_vir_frac):
        _lay_v = _spec.setdefault("layout", {})
        _lay_v.pop("config", None)
        _lay_v["dragmode"] = False
    vir_reads_json = json.dumps(spec_vir_reads, cls=PlotlyJSONEncoder, allow_nan=False)
    vir_frac_json = json.dumps(spec_vir_frac, cls=PlotlyJSONEncoder, allow_nan=False)

    page = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8"/>
  <meta name="viewport" content="width=device-width, initial-scale=1"/>
  <title>Metatropics — Read Distribution</title>
  <script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
  <style>
    :root {{
      --bg: #f1f5f9;
      --card: #ffffff;
      --border: #e2e8f0;
      --muted: #475569;
      --text: #0f172a;
      --accent: #4f46e5;
    }}
    body {{
      margin: 0;
      font-family: system-ui, Segoe UI, Roboto, Helvetica, Arial, sans-serif;
      background: linear-gradient(165deg, #e8eef7 0%, var(--bg) 42%, #f8fafc 100%);
      color: var(--text);
      min-height: 100vh;
    }}
    .wrap {{
      max-width: 1280px;
      margin: 0 auto;
      padding: 28px 20px 40px;
    }}
    .card {{
      background: var(--card);
      border: 1px solid var(--border);
      border-radius: 14px;
      box-shadow: 0 10px 40px rgba(15, 23, 42, 0.08);
      padding: 16px 16px 10px;
    }}
    .chart-title {{
      margin: 0 0 14px;
      font-size: 1.45rem;
      font-weight: 700;
      letter-spacing: -0.02em;
      color: var(--text);
    }}
    .toolbar {{
      display: flex;
      flex-wrap: wrap;
      gap: 10px;
      align-items: center;
      margin: 0 0 12px;
    }}
    .seg {{
      display: inline-flex;
      border: 1px solid var(--border);
      border-radius: 999px;
      overflow: hidden;
      background: #f8fafc;
    }}
    .seg button {{
      border: 0;
      background: transparent;
      color: var(--muted);
      padding: 8px 14px;
      font-size: 0.9rem;
      cursor: pointer;
    }}
    .seg button.active {{
      background: rgba(79, 70, 229, 0.12);
      color: var(--accent);
      font-weight: 600;
    }}
    .dl {{
      border: 1px solid var(--border);
      background: #fff;
      color: var(--text);
      border-radius: 999px;
      padding: 7px 12px;
      font-size: 0.85rem;
      cursor: pointer;
    }}
    .dl:hover {{
      border-color: #cbd5e1;
      background: #f8fafc;
    }}
    #plot-pct, #plot-cnt, #plot-vir {{
      width: 100%;
      min-height: 420px;
    }}
    /* Nudge Plotly mode bar slightly up (zoom / autoscale); PNG export is the toolbar button only. */
    #plot-pct .modebar,
    #plot-cnt .modebar,
    #plot-vir .modebar {{
      top: 2px !important;
      right: 8px !important;
    }}
  </style>
</head>
<body>
  <div class="wrap">
    <div class="card">
      <h1 class="chart-title">Read Distribution</h1>
      <div class="toolbar">
        <div class="seg" role="group" aria-label="Chart scale">
          <button type="button" id="btn-pct" class="active" onclick="showPct()">Percent</button>
          <button type="button" id="btn-cnt" onclick="showCnt()">Counts</button>
          <button type="button" id="btn-vir" onclick="showVir()">Viral detail</button>
        </div>
        <button type="button" class="dl" onclick="dlActive('png')">PNG</button>
        <button type="button" class="dl" onclick="dlActive('svg')">SVG</button>
      </div>
      <div id="vir-toolbar" class="toolbar" style="display:none">
        <div class="seg" role="group" aria-label="Viral detail view">
          <button type="button" id="btn-vir-reads" class="active" onclick="showVirReads()">Reads</button>
          <button type="button" id="btn-vir-pct" onclick="showVirFrac()">100%</button>
        </div>
      </div>
      <div id="plot-pct"></div>
      <div id="plot-cnt" style="display:none"></div>
      <div id="plot-vir" style="display:none"></div>
    </div>
  </div>
  <script>
    const specPct = {pct_json};
    const specCnt = {cnt_json};
    const specVirReads = {vir_reads_json};
    const specVirFrac = {vir_frac_json};
    const cfg = {{
      responsive: true,
      displayModeBar: true,
      displaylogo: false,
      scrollZoom: false,
      // PNG uses Plotly.downloadImage from the left "PNG" button (scale 3). No camera on the mode bar.
      modeBarButtons: [[
        'zoomIn2d',
        'zoomOut2d',
        'autoScale2d',
        'resetScale2d'
      ]]
    }};
    function cleanLayout(layout) {{
      if (!layout) return;
      delete layout.config;
      layout.dragmode = false;
    }}
    function mount(id, spec) {{
      cleanLayout(spec.layout);
      Plotly.newPlot(id, spec.data, spec.layout, cfg);
    }}
    cleanLayout(specVirReads.layout);
    cleanLayout(specVirFrac.layout);
    mount('plot-pct', specPct);
    mount('plot-cnt', specCnt);
    mount('plot-vir', specVirReads);
    function activePlotId() {{
      if (document.getElementById('plot-vir').style.display !== 'none') return 'plot-vir';
      if (document.getElementById('plot-cnt').style.display !== 'none') return 'plot-cnt';
      return 'plot-pct';
    }}
    function dlActive(fmt) {{
      const id = activePlotId();
      let base = 'read_distribution_percent';
      if (id === 'plot-cnt') base = 'read_distribution_counts';
      if (id === 'plot-vir') {{
        const pct = document.getElementById('btn-vir-pct').classList.contains('active');
        base = pct ? 'read_distribution_viral_100pct' : 'read_distribution_viral_reads';
      }}
      const opts = {{ format: fmt, filename: base }};
      if (fmt === 'png') opts.scale = 3;
      Plotly.downloadImage(document.getElementById(id), opts);
    }}
    function showPct() {{
      document.getElementById('plot-pct').style.display = 'block';
      document.getElementById('plot-cnt').style.display = 'none';
      document.getElementById('plot-vir').style.display = 'none';
      document.getElementById('vir-toolbar').style.display = 'none';
      document.getElementById('btn-pct').classList.add('active');
      document.getElementById('btn-cnt').classList.remove('active');
      document.getElementById('btn-vir').classList.remove('active');
      Plotly.Plots.resize(document.getElementById('plot-pct'));
    }}
    function showCnt() {{
      document.getElementById('plot-pct').style.display = 'none';
      document.getElementById('plot-cnt').style.display = 'block';
      document.getElementById('plot-vir').style.display = 'none';
      document.getElementById('vir-toolbar').style.display = 'none';
      document.getElementById('btn-cnt').classList.add('active');
      document.getElementById('btn-pct').classList.remove('active');
      document.getElementById('btn-vir').classList.remove('active');
      Plotly.Plots.resize(document.getElementById('plot-cnt'));
    }}
    function showVirReads() {{
      document.getElementById('btn-vir-reads').classList.add('active');
      document.getElementById('btn-vir-pct').classList.remove('active');
      Plotly.react(document.getElementById('plot-vir'), specVirReads.data, specVirReads.layout, cfg);
    }}
    function showVirFrac() {{
      document.getElementById('btn-vir-pct').classList.add('active');
      document.getElementById('btn-vir-reads').classList.remove('active');
      Plotly.react(document.getElementById('plot-vir'), specVirFrac.data, specVirFrac.layout, cfg);
    }}
    function showVir() {{
      document.getElementById('plot-pct').style.display = 'none';
      document.getElementById('plot-cnt').style.display = 'none';
      document.getElementById('plot-vir').style.display = 'block';
      document.getElementById('vir-toolbar').style.display = 'flex';
      document.getElementById('btn-vir').classList.add('active');
      document.getElementById('btn-pct').classList.remove('active');
      document.getElementById('btn-cnt').classList.remove('active');
      showVirReads();
      Plotly.Plots.resize(document.getElementById('plot-vir'));
    }}
    window.addEventListener('resize', function() {{
      try {{
        Plotly.Plots.resize(document.getElementById('plot-pct'));
        Plotly.Plots.resize(document.getElementById('plot-cnt'));
        Plotly.Plots.resize(document.getElementById('plot-vir'));
      }} catch (e) {{}}
    }});
  </script>
</body>
</html>
"""
    html_out.write_text(page, encoding="utf-8")


def render_pdf_matplotlib(pdf_out: Path, csv_path: Path) -> None:
    import pandas as pd  # type: ignore
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt  # type: ignore

    df = pd.read_csv(csv_path)
    pct_cols = [f"{c}_pct" for c in CATEGORY_ORDER if f"{c}_pct" in df.columns]
    if not pct_cols:
        raise RuntimeError("No *_pct columns found in read_counts.csv")

    if "viral_pct" in df.columns:
        df = df.sort_values(by=["viral_pct", "sample"], ascending=[False, True]).reset_index(drop=True)

    cat_order = pct_cols
    labels = [CATEGORY_LABELS.get(c.replace("_pct", ""), c.replace("_pct", "")) for c in cat_order]
    cols = [CATEGORY_COLORS.get(c.replace("_pct", ""), "#94a3b8") for c in cat_order]

    try:
        plt.style.use("seaborn-v0_8-white")
    except Exception:
        try:
            plt.style.use("seaborn-white")
        except Exception:
            pass

    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica", "Liberation Sans"],
            "font.size": 9,
            "axes.titlesize": 11,
            "axes.labelsize": 10,
            "legend.fontsize": 8,
            "axes.titleweight": "bold",
        }
    )

    n = len(df)
    row_h = 0.32
    fig_h = max(4.2, min(22.0, 1.15 + row_h * n))
    fig_w = 7.2  # single-column width (inches)

    fig, ax = plt.subplots(figsize=(fig_w, fig_h), dpi=300)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("#fafafa")

    y = list(range(n))
    left = [0.0] * n
    for c, label, color in zip(cat_order, labels, cols):
        vals = df[c].fillna(0).astype(float).tolist()
        ax.barh(
            y,
            vals,
            left=left,
            label=label,
            color=color,
            edgecolor="white",
            linewidth=0.6,
            height=0.72,
        )
        left = [l + v for l, v in zip(left, vals)]

    ax.set_yticks(y)
    ax.set_yticklabels(df["sample"].astype(str).tolist(), fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("Percentage (%)")
    ax.set_title("Read Distribution")
    ax.set_xlim(0, 100)
    ax.xaxis.set_major_locator(plt.MultipleLocator(10))
    ax.grid(axis="x", linestyle="-", linewidth=0.35, color="#d4d4d8", alpha=0.9)
    ax.grid(axis="y", visible=False)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.spines["left"].set_color("#cbd5e1")
    ax.spines["bottom"].set_color("#cbd5e1")

    # One legend row with one column per category so entries span the figure width (no empty corner).
    ncol = len(labels) if len(labels) <= 6 else min(3, len(labels))
    leg = ax.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, -0.16),
        ncol=ncol,
        frameon=True,
        fancybox=False,
        edgecolor="#e2e8f0",
        columnspacing=1.35 if ncol == len(labels) else 1.05,
        handlelength=1.0,
        handletextpad=0.5,
        borderaxespad=1.35,
        borderpad=0.38,
        labelspacing=0.45,
        title="Category",
        title_fontsize=9,
        fontsize=8,
    )

    n_rows = max(1, (len(labels) + ncol - 1) // ncol)
    # Extra room for legend title blank line above category swatches
    bottom = min(0.52, 0.18 + 0.065 * n_rows)
    fig.subplots_adjust(left=0.22, right=0.98, top=0.92, bottom=bottom)
    extra = [leg] if leg is not None else []
    fig.savefig(pdf_out, format="pdf", bbox_extra_artists=extra, bbox_inches="tight", pad_inches=0.16, dpi=300)
    plt.close(fig)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", required=True, help="Pipeline outdir (contains Classification/virasign, Reads/* etc)")
    ap.add_argument("--host-status", default="not_used", help="not_used|human_only|other_only|both")
    ap.add_argument("--workdir", default=".", help="Directory containing read_count/ staging folder")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    workdir = Path(args.workdir)
    read_count_dir = workdir / "read_count"
    read_count_dir.mkdir(parents=True, exist_ok=True)

    rows, include_human, include_other = compute_rows(read_count_dir, outdir, args.host_status)
    out_csv = read_count_dir / "read_counts.csv"
    write_csv(rows, include_human, include_other, out_csv)

    html_out = read_count_dir / "read_distribution.html"
    pdf_out = read_count_dir / "read_distribution.pdf"
    render_plotly_html(html_out, out_csv, outdir)
    render_pdf_matplotlib(pdf_out, out_csv)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

