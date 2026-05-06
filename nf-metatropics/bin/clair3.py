#!/usr/bin/env python3
"""
Uniform Clair3 post-processing for Metatropics.

Purpose:
- Re-count ref/alt support from the per-virus BAM (SNPs + simple left-anchored indels).
- Apply uniform filters (QUAL when present, DP, ALT reads, optional strand-bias p-value).
- Apply MAPQ filter during recount (min_mq).
- Classify variants as major or minor by VAF thresholds (stored as INFO/TIER=major|minor in VCF).
- Emit:
  - <prefix>.variants.filtered.vcf
  - <prefix>.variants.unfiltered.vcf  (includes INFO/UFR and INFO/FIN)
  - <prefix>.variants.html

This script is adapted from the original Metatropics Medaka post-processor, but uses Clair3 VCF as input.

Dependencies: pysam.
"""

from __future__ import annotations

import argparse
import gzip
import html
import math
import re
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pysam


@dataclass(frozen=True)
class VcfCall:
    chrom: str
    pos: int
    ref: str
    alt: str
    qual: Optional[float]
    filt: str
    info: Dict[str, str]
    format_keys: Optional[str]
    sample_cell: Optional[str]


def _safe_float(x: Optional[str]) -> Optional[float]:
    if x is None:
        return None
    try:
        if x == ".":
            return None
        return float(x)
    except Exception:
        return None


def _parse_info(info: str) -> Dict[str, str]:
    out: Dict[str, str] = {}
    if not info or info == ".":
        return out
    for kv in info.split(";"):
        if not kv:
            continue
        if "=" in kv:
            k, v = kv.split("=", 1)
            out[k] = v
        else:
            out[kv] = "true"
    return out


def _read_vcf_calls(path: Path) -> Tuple[List[str], List[VcfCall]]:
    header_lines: List[str] = []
    calls: List[VcfCall] = []
    opener = gzip.open if path.name.endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line:
                continue
            if line.startswith("#"):
                header_lines.append(line.rstrip("\n"))
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            chrom, pos_s, _id, ref, alt, qual_s, flt, info = parts[:8]
            fmt = parts[8] if len(parts) > 8 else None
            smp = parts[9] if len(parts) > 9 else None
            info_dict = _parse_info(info)
            for a in alt.split(","):
                calls.append(
                    VcfCall(
                        chrom=chrom,
                        pos=int(pos_s),
                        ref=ref,
                        alt=a,
                        qual=_safe_float(qual_s),
                        filt=flt,
                        info=info_dict,
                        format_keys=fmt,
                        sample_cell=smp,
                    )
                )
    return header_lines, calls


def _is_snp(ref: str, alt: str) -> bool:
    return len(ref) == 1 and len(alt) == 1 and ref != alt and ref != "-" and alt != "-"


def _is_indel(ref: str, alt: str) -> bool:
    if ref == "-" or alt == "-":
        return True
    return len(ref) != len(alt)


def _count_alleles_at_variant(
    bam: pysam.AlignmentFile,
    fasta: pysam.FastaFile,
    chrom: str,
    pos1: int,
    ref: str,
    alt: str,
    min_bq: int,
    min_mq: int,
) -> Dict[str, int]:
    pos0 = pos1 - 1

    is_snp = _is_snp(ref, alt)
    is_indel = _is_indel(ref, alt)

    ref_f = ref_r = alt_f = alt_r = 0
    dp = 0

    for pileup_col in bam.pileup(
        chrom,
        pos0,
        pos0 + 1,
        truncate=True,
        stepper="nofilter",
        min_base_quality=0,
        ignore_orphans=False,
        max_depth=500000,
    ):
        if pileup_col.reference_pos != pos0:
            continue
        for pr in pileup_col.pileups:
            aln = pr.alignment
            if aln.is_unmapped or aln.is_secondary or aln.is_supplementary:
                continue
            if aln.mapping_quality < min_mq:
                continue

            is_rev = aln.is_reverse

            if is_snp:
                qpos = pr.query_position
                if qpos is None:
                    continue
                bq = aln.query_qualities[qpos] if aln.query_qualities is not None else 60
                if bq < min_bq:
                    continue
                base = aln.query_sequence[qpos]
                if base == ref:
                    dp += 1
                    if is_rev:
                        ref_r += 1
                    else:
                        ref_f += 1
                elif base == alt:
                    dp += 1
                    if is_rev:
                        alt_r += 1
                    else:
                        alt_f += 1
                else:
                    continue

            elif is_indel:
                qpos = pr.query_position
                if qpos is None:
                    continue
                bq = aln.query_qualities[qpos] if aln.query_qualities is not None else 60
                if bq < min_bq:
                    continue

                exp_indel_len = len(alt) - len(ref)
                obs_indel_len = pr.indel

                base = aln.query_sequence[qpos]
                if obs_indel_len == 0:
                    if base == ref[0]:
                        dp += 1
                        if is_rev:
                            ref_r += 1
                        else:
                            ref_f += 1
                    continue

                if obs_indel_len != exp_indel_len:
                    continue

                if exp_indel_len > 0:
                    ins_seq = aln.query_sequence[qpos + 1 : qpos + 1 + exp_indel_len]
                    exp_seq = alt[len(ref) :]
                    if base != ref[0]:
                        continue
                    if ins_seq != exp_seq:
                        continue
                    dp += 1
                    if is_rev:
                        alt_r += 1
                    else:
                        alt_f += 1
                    continue

                if exp_indel_len < 0:
                    del_len = -exp_indel_len
                    exp_deleted = ref[len(alt) :]
                    ref_deleted = fasta.fetch(chrom, pos0 + 1, pos0 + 1 + del_len)
                    if base != ref[0]:
                        continue
                    if ref_deleted != exp_deleted:
                        continue
                    dp += 1
                    if is_rev:
                        alt_r += 1
                    else:
                        alt_f += 1
                    continue

    return {"DP": dp, "REF_F": ref_f, "REF_R": ref_r, "ALT_F": alt_f, "ALT_R": alt_r}


def _log_factorial(n: int) -> float:
    if n <= 1:
        return 0.0
    return math.lgamma(n + 1)


def _log_binom_pmf_half(n: int, k: int) -> float:
    return _log_factorial(n) - _log_factorial(k) - _log_factorial(n - k) - n * math.log(2.0)


def _binom_twosided_strand_alt(alt_f: int, alt_r: int) -> Optional[float]:
    n = alt_f + alt_r
    if n <= 0 or alt_f < 0 or alt_r < 0:
        return None
    if n > 8000:
        z = (alt_f - n * 0.5) / math.sqrt(n * 0.25 + 1e-12)
        return float(min(1.0, max(0.0, math.erfc(abs(z) / math.sqrt(2.0)))))
    log_obs = _log_binom_pmf_half(n, alt_f)
    p = 0.0
    for x in range(n + 1):
        lp = _log_binom_pmf_half(n, x)
        if lp <= log_obs + 1e-14:
            p += math.exp(lp)
    return float(min(1.0, max(0.0, p)))


def _strand_bias_pvalue(ref_f: int, ref_r: int, alt_f: int, alt_r: int) -> Optional[float]:
    _ = (ref_f, ref_r)
    return _binom_twosided_strand_alt(alt_f, alt_r)


def _apobec_context(
    fasta: pysam.FastaFile,
    chrom: str,
    pos1: int,
    ref: str,
    alt: str,
) -> Tuple[Optional[str], Optional[str], bool, bool]:
    if not _is_snp(ref, alt):
        return None, None, False, False

    pos0 = pos1 - 1
    ctx = fasta.fetch(chrom, pos0 - 1, pos0 + 2).upper()
    if len(ctx) != 3 or "N" in ctx:
        r = ref.upper()
        a = alt.upper()
        is_potential = (sorted([r, a]) == ["C", "T"]) or (sorted([r, a]) == ["A", "G"])
        return ctx if ctx else None, None, False, bool(is_potential)

    r = ref.upper()
    a = alt.upper()
    is_ct = sorted([r, a]) == ["C", "T"]
    is_ag = sorted([r, a]) == ["A", "G"]
    is_potential = is_ct or is_ag

    is_apobec = False
    apobec_tag: Optional[str] = None

    if is_ct:
        if r == "C" and a == "T" and ctx[0] == "T" and ctx[1] == "C":
            is_apobec, apobec_tag = True, "TC_TT"
        elif r == "T" and a == "C" and ctx[0] == "T" and ctx[1] == "T":
            is_apobec, apobec_tag = True, "TT_TC"

    if is_ag:
        if r == "G" and a == "A" and ctx[1] == "G" and ctx[2] == "A":
            is_apobec, apobec_tag = True, "GA_AA"
        elif r == "A" and a == "G" and ctx[1] == "A" and ctx[2] == "A":
            is_apobec, apobec_tag = True, "AA_GA"

    if is_apobec:
        return ctx, apobec_tag, True, is_potential
    if is_ct:
        return ctx, "C<->T_other", False, is_potential
    if is_ag:
        return ctx, "A<->G_other", False, is_potential
    return ctx, None, False, False


_INFO_LINE_RE = re.compile(r"^##INFO=<ID=([^,>]+)")


def _merge_header_lines(base_header_lines: List[str], extra_info_defs: List[str]) -> List[str]:
    existing_ids = set()
    out: List[str] = []
    for line in base_header_lines:
        if line.startswith("##INFO=<ID="):
            m = _INFO_LINE_RE.match(line)
            if m:
                existing_ids.add(m.group(1))
        out.append(line)

    insert_at = max(0, len(out) - 1)
    for d in extra_info_defs:
        m = _INFO_LINE_RE.match(d)
        if not m:
            continue
        iid = m.group(1)
        if iid in existing_ids:
            continue
        out.insert(insert_at, d)
        insert_at += 1
        existing_ids.add(iid)
    return out


def _fmt_info(kv: Dict[str, object]) -> str:
    parts: List[str] = []
    for k in sorted(kv.keys()):
        v = kv[k]
        if v is True:
            parts.append(k)
        elif v is False or v is None:
            continue
        else:
            parts.append(f"{k}={v}")
    return ";".join(parts)


def _pct(x: object) -> str:
    try:
        return f"{float(x) * 100:.0f}%"
    except (TypeError, ValueError):
        return html.escape(str(x))


def _build_settings_inner_html(settings: dict) -> str:
    uf = settings.get("uniform_filters") or {}
    min_sb = float(uf.get("min_sb_pvalue", 0) or 0)

    major_v = uf.get("major_vaf")
    mmin, mmax = uf.get("minor_vaf_min"), uf.get("minor_vaf_max")
    minor_band = f"{_pct(mmin)}-{_pct(mmax)}"

    acc = str(settings.get("virus") or "").strip()
    ref_box: Optional[Tuple[str, str]] = ("Reference", html.escape(acc)) if acc else None

    boxes: List[Tuple[str, str]] = [
        ("Min quality", html.escape(str(uf.get("min_qual", "")))),
        ("Min depth (DP)", html.escape(str(uf.get("min_dp", "")))),
        ("Min ALT reads", html.escape(str(uf.get("min_alt_reads", "")))),
        ("Min MAPQ", html.escape(str(uf.get("min_mq", "")))),
        ("Major", html.escape(f"VAF ≥ {_pct(major_v)}")),
        ("Minor", minor_band),
    ]
    # Only show a Fisher/SB threshold box if enabled.
    # If disabled, use that slot for the Reference (so we never show "off").
    if min_sb > 0:
        boxes.append(("SB Fisher", f"&gt; {min_sb:g}"))
        if ref_box:
            boxes.append(ref_box)
    else:
        if ref_box:
            boxes.append(ref_box)

    bits = []
    for label, val in boxes:
        bits.append(
            f'<div class="stat-box"><div class="stat-label">{html.escape(label)}</div>'
            f'<div class="stat-value">{val}</div></div>'
        )
    return f'<div class="settings-row"><div class="stat-grid">{"".join(bits)}</div></div>'


def _tier_cell_display(tier: str) -> str:
    if tier == "major":
        return "Major"
    if tier == "minor":
        return "Minor"
    return tier


def _build_variant_table_fragment(
    rows: List[dict],
    section_title: str,
    *,
    include_ctx: bool,
    include_apobec3: bool,
) -> str:
    cols: List[Tuple[str, str]] = [
        ("POS", "POS"),
        ("REF", "REF"),
        ("ALT", "ALT"),
        ("TIER", "Variants"),
        ("QUAL", "QUAL"),
        ("DP", "DP"),
        ("AD", "AD"),
        ("VAF", "VAF"),
        ("SB", "SB"),
    ]
    if include_ctx:
        cols.append(("CTX", "CTX"))
    if include_apobec3:
        cols.append(("APOBEC3", "APOBEC3"))
    thead = '<th class="nr">Nr</th>' + "".join(f"<th>{html.escape(c[1])}</th>" for c in cols)
    if not rows:
        return f"""        <h3 class="section-title">{html.escape(section_title)}</h3>
        <p class="settings-prose">None in this call set.</p>
"""
    tbody_lines: List[str] = []
    for nr, r in enumerate(rows, start=1):
        tds = [f'<td class="nr">{nr}</td>']
        for key, _label in cols:
            v = r.get(key, "")
            if key == "TIER":
                v = _tier_cell_display(str(v))
            elif key == "SB" and (v is None or str(v).strip() == ""):
                v = "NA"
            if v is None:
                v = ""
            tds.append(f"<td>{html.escape(str(v))}</td>")
        tbody_lines.append("<tr>" + "".join(tds) + "</tr>")
    tbody = "\n".join(tbody_lines)
    return f"""        <h3 class="section-title">{html.escape(section_title)}</h3>
        <div class="table-wrap">
          <table>
            <thead><tr>{thead}</tr></thead>
            <tbody>
{tbody}
            </tbody>
          </table>
        </div>
"""


def _build_html_report(*, sample: str, virus: str, settings: dict, rows: List[dict], out_path: Path) -> None:
    snp_rows = [r for r in rows if r.get("TYPE") == "SNP"]
    indel_rows = [r for r in rows if r.get("TYPE") != "SNP"]
    snp_rows.sort(key=lambda r: (int(r["POS"]), str(r["REF"]), str(r["ALT"])))
    indel_rows.sort(key=lambda r: (int(r["POS"]), str(r["REF"]), str(r["ALT"])))

    n_snp, n_indel = len(snp_rows), len(indel_rows)

    summary_line = (
        f'<p class="variant-count"><strong>{n_snp}</strong> SNPs · <strong>{n_indel}</strong> indels.</p>'
        if (n_snp or n_indel)
        else '<p class="variant-count">No variant sites passed filters and major/minor classification.</p>'
    )

    qual_note = (
        '<p class="qual-hint"><strong>QUAL</strong> is the variant-level confidence from the caller (Phred-like: '
        "higher = more confident). It is not FASTQ base quality.</p>"
    )

    settings_boxes = _build_settings_inner_html(settings)
    snp_block = _build_variant_table_fragment(snp_rows, "SNPs", include_ctx=True, include_apobec3=True)
    indel_block = _build_variant_table_fragment(indel_rows, "Indels", include_ctx=False, include_apobec3=False)

    doc = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8"/>
  <meta name="viewport" content="width=device-width, initial-scale=1"/>
  <title>Clair3 variants — {html.escape(sample)} / {html.escape(virus)}</title>
  <style>
    :root {{
      --bg: #0b1020;
      --panel: #121a33;
      --text: #e8ecff;
      --muted: #a9b4e6;
      --border: rgba(255,255,255,0.10);
    }}
    body {{
      margin: 0;
      font-family: ui-sans-serif, system-ui, -apple-system, Segoe UI, Roboto, Arial, sans-serif;
      background: radial-gradient(1200px 800px at 20% 0%, #1b2a66 0%, var(--bg) 55%);
      color: var(--text);
    }}
    .wrap {{ max-width: 1200px; margin: 0 auto; padding: 28px 20px 40px; }}
    h1 {{ font-size: 22px; margin: 0 0 20px; letter-spacing: 0.2px; }}
    .grid {{ display: grid; grid-template-columns: 1fr; gap: 14px; }}
    .card {{
      background: linear-gradient(180deg, rgba(255,255,255,0.06), rgba(255,255,255,0.03));
      border: 1px solid var(--border);
      border-radius: 14px;
      padding: 14px 14px;
      box-shadow: 0 10px 30px rgba(0,0,0,0.25);
    }}
    .pill {{
      display: inline-flex; gap: 8px; align-items: center;
      padding: 6px 10px; border-radius: 999px;
      border: 1px solid var(--border);
      background: rgba(0,0,0,0.18);
      color: var(--muted);
      font-size: 12px;
    }}
    .settings-prose {{
      margin: 0;
      padding: 2px 0 0;
      color: #dbe4ff;
      font-size: 14px;
      line-height: 1.55;
    }}
    .settings-row {{
      display: flex;
      flex-wrap: wrap;
      align-items: flex-start;
      justify-content: space-between;
      gap: 14px 20px;
      margin-top: 4px;
    }}
    .stat-grid {{
      display: grid;
      grid-template-columns: repeat(auto-fill, minmax(148px, 1fr));
      gap: 10px;
      flex: 1 1 280px;
      min-width: 0;
      width: 100%;
    }}
    .stat-box {{
      border: 1px solid var(--border);
      border-radius: 10px;
      padding: 10px 12px;
      background: rgba(0,0,0,0.22);
    }}
    .stat-label {{
      font-size: 11px;
      text-transform: uppercase;
      letter-spacing: 0.06em;
      color: var(--muted);
      margin-bottom: 6px;
    }}
    .stat-value {{
      font-size: 17px;
      font-weight: 650;
      color: #dbe4ff;
      line-height: 1.25;
    }}
    .section-title {{
      font-size: 15px;
      font-weight: 650;
      margin: 18px 0 10px;
      color: #e8ecff;
    }}
    .qual-hint {{
      margin: 14px 0 0;
      font-size: 12px;
      color: var(--muted);
      line-height: 1.45;
    }}
    .variant-count {{
      margin: 0 0 12px;
      font-size: 13px;
      color: var(--muted);
      line-height: 1.45;
    }}
    th.nr, td.nr {{
      width: 3.2rem;
      text-align: right;
      font-variant-numeric: tabular-nums;
      color: #a9b4e6;
    }}
    .table-wrap {{
      overflow: auto;
      border-radius: 12px;
      border: 1px solid var(--border);
    }}
    table {{ border-collapse: separate; border-spacing: 0; width: 100%; min-width: 980px; }}
    th, td {{ padding: 10px 10px; border-bottom: 1px solid var(--border); vertical-align: top; }}
    th {{
      position: sticky; top: 0;
      background: rgba(18, 26, 51, 0.92);
      backdrop-filter: blur(8px);
      text-align: left;
      font-size: 12px;
      letter-spacing: 0.6px;
      text-transform: uppercase;
      color: #cfd8ff;
      border-bottom: 1px solid rgba(255,255,255,0.14);
    }}
    tr:hover td {{ background: rgba(122, 162, 255, 0.06); }}
    td {{ font-size: 13px; color: #eef1ff; }}
  </style>
</head>
<body>
  <div class="wrap">
    <div class="pill">Metatropics · Clair3 variant calling</div>
    <h1>Variants</h1>
    <div class="grid">
      <div class="card">
        <div style="font-weight:650; margin-bottom:8px;">Settings</div>
        {settings_boxes}
      </div>
      <div class="card">
        {summary_line}
{snp_block}
{indel_block}
        {qual_note}
      </div>
    </div>
  </div>
</body>
</html>
"""
    out_path.write_text(doc, encoding="utf-8")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--virus", required=True)
    ap.add_argument("--bam", required=True, type=Path)
    ap.add_argument("--ref-fasta", required=True, type=Path)
    ap.add_argument("--clair3-vcf", required=True, type=Path)
    ap.add_argument("--out-prefix", required=True, help="Prefix for outputs in cwd")

    ap.add_argument("--min-qual", type=float, required=True)
    ap.add_argument("--min-bq", type=int, required=True)
    ap.add_argument("--min-mq", type=int, required=True)
    ap.add_argument("--min-dp", type=int, required=True)
    ap.add_argument("--min-alt-reads", type=int, required=True)
    ap.add_argument("--major-vaf", type=float, required=True)
    ap.add_argument("--minor-vaf-min", type=float, required=True)
    ap.add_argument("--minor-vaf-max", type=float, required=True)
    ap.add_argument("--min-sb-pvalue", type=float, required=True)
    ap.add_argument(
        "--sb-min-alt-strand",
        type=int,
        default=0,
        help=(
            "Skip the SB (strand-bias) filter for non-APOBEC3 sites when BOTH ALT_F and ALT_R "
            "are >= this value. Set to 0 to always apply SB filtering when enabled."
        ),
    )
    args = ap.parse_args()

    prefix = str(args.out_prefix)

    fasta = pysam.FastaFile(str(args.ref_fasta))
    bam = pysam.AlignmentFile(str(args.bam), "rb")

    contigs = list(fasta.references)
    if len(contigs) != 1:
        raise SystemExit(f"Expected exactly 1 contig in {args.ref_fasta}, got {contigs}")
    contig = contigs[0]
    contig_len = int(fasta.get_reference_length(contig))

    header_lines, calls = _read_vcf_calls(args.clair3_vcf)
    if not header_lines or not any(h.startswith("#CHROM") for h in header_lines):
        raise SystemExit(f"Missing VCF header in {args.clair3_vcf}")

    extra_info_defs = [
        '##INFO=<ID=TIER,Number=1,Type=String,Description="Allele fraction class (values major|minor; major or minor variant by VAF band)">',
        '##INFO=<ID=TYPE,Number=1,Type=String,Description="Variant type: SNP|INDEL">',
        '##INFO=<ID=DP,Number=1,Type=Integer,Description="Biallelic ref+alt supporting read depth (uniform recount, MAPQ-filtered)">',
        '##INFO=<ID=AD,Number=2,Type=Integer,Description="Allelic depths (REF,ALT) from uniform recount (MAPQ-filtered)">',
        '##INFO=<ID=VAF,Number=1,Type=Float,Description="ALT/(REF+ALT) from uniform recount (MAPQ-filtered)">',
        '##INFO=<ID=ALT_F,Number=1,Type=Integer,Description="ALT supporting reads on forward strand (uniform recount, MAPQ-filtered)">',
        '##INFO=<ID=ALT_R,Number=1,Type=Integer,Description="ALT supporting reads on reverse strand (uniform recount, MAPQ-filtered)">',
        '##INFO=<ID=REF_F,Number=1,Type=Integer,Description="REF supporting reads on forward strand (uniform recount, MAPQ-filtered)">',
        '##INFO=<ID=REF_R,Number=1,Type=Integer,Description="REF supporting reads on reverse strand (uniform recount, MAPQ-filtered)">',
        '##INFO=<ID=SB,Number=1,Type=Float,Description="ALT strand-bias p-value (two-sided exact): forward vs reverse among ALT-supporting reads (Binomial n=ALT_F+ALT_R, p=0.5); REF reads not used">',
        '##INFO=<ID=CTX,Number=1,Type=String,Description="Reference 3-mer context centered on POS (plus strand)">',
        '##INFO=<ID=APOBEC3,Number=1,Type=String,Description="APOBEC3 dinucleotide motif on ref: C<->T with 5-prime T (TC/TT) or G<->A with 3-prime A (GA/AA); Yes or No">',
    ]
    extra_info_unf_only = [
        '##INFO=<ID=UFR,Number=1,Type=String,Description="Uniform filter status: PASS or comma-separated failure codes">',
        '##INFO=<ID=FIN,Number=1,Type=Integer,Description="1 if row is emitted in variants.filtered VCF else 0">',
    ]

    merged_header = _merge_header_lines(header_lines, extra_info_defs)
    merged_header_unf = _merge_header_lines(header_lines, extra_info_defs + extra_info_unf_only)

    shutil.copy2(args.clair3_vcf, Path(f"{prefix}.clair3.in.vcf.gz"))

    rows_out: List[dict] = []
    vcf_body_lines: List[str] = []
    unf_body_lines: List[str] = []

    for c in calls:
        if c.chrom != contig:
            continue

        counts = _count_alleles_at_variant(
            bam,
            fasta,
            c.chrom,
            c.pos,
            c.ref,
            c.alt,
            min_bq=int(args.min_bq),
            min_mq=int(args.min_mq),
        )
        dp = int(counts["DP"])
        ad_ref = int(counts["REF_F"] + counts["REF_R"])
        ad_alt = int(counts["ALT_F"] + counts["ALT_R"])
        alt_f = int(counts["ALT_F"])
        alt_r = int(counts["ALT_R"])
        if dp <= 0:
            continue

        vaf = float(ad_alt / dp) if dp else 0.0
        sb = _strand_bias_pvalue(int(counts["REF_F"]), int(counts["REF_R"]), alt_f, alt_r)
        ctx, _, is_apobec, _ = _apobec_context(fasta, c.chrom, c.pos, c.ref, c.alt)
        vtype = "SNP" if _is_snp(c.ref, c.alt) else "INDEL"

        fail_chain: List[str] = []
        if dp < int(args.min_dp):
            fail_chain.append("min_dp")
        if ad_alt < int(args.min_alt_reads):
            fail_chain.append("min_alt_reads")
        if c.qual is not None and float(c.qual) < float(args.min_qual):
            fail_chain.append("min_qual")
        if (
            float(args.min_sb_pvalue) > 0.0
            and sb is not None
            and float(sb) < float(args.min_sb_pvalue)
            and not is_apobec
            and not (
                int(args.sb_min_alt_strand) > 0
                and alt_f >= int(args.sb_min_alt_strand)
                and alt_r >= int(args.sb_min_alt_strand)
            )
        ):
            fail_chain.append("strand_bias")

        passes_uniform = len(fail_chain) == 0
        tier_label = "filtered"
        in_tiered_vcf = False
        if passes_uniform:
            if vaf >= float(args.major_vaf):
                tier_label = "major"
                in_tiered_vcf = True
            elif float(args.minor_vaf_min) <= vaf < float(args.minor_vaf_max):
                tier_label = "minor"
                in_tiered_vcf = True
            else:
                tier_label = "other_vaf"

        ufr = ",".join(fail_chain) if fail_chain else "PASS"
        qual_out = "." if c.qual is None else f"{c.qual:.3g}"
        fmt = c.format_keys if c.format_keys else "GT"
        smp = c.sample_cell if c.sample_cell else "1/1"

        unf_info: Dict[str, object] = dict(c.info)
        unf_info.update(
            {
                "TIER": tier_label,
                "TYPE": vtype,
                "DP": dp,
                "AD": f"{ad_ref},{ad_alt}",
                "VAF": f"{vaf:.5f}",
                "ALT_F": alt_f,
                "ALT_R": alt_r,
                "REF_F": int(counts["REF_F"]),
                "REF_R": int(counts["REF_R"]),
                "SB": f"{sb:.3g}" if sb is not None else ".",
                "CTX": ctx if ctx else ".",
                "APOBEC3": "Yes" if is_apobec else "No",
                "UFR": ufr,
                "FIN": 1 if in_tiered_vcf else 0,
            }
        )
        unf_body_lines.append(
            "\t".join(
                [
                    c.chrom,
                    str(c.pos),
                    ".",
                    c.ref,
                    c.alt,
                    qual_out,
                    "PASS",
                    _fmt_info(unf_info),
                    fmt,
                    smp,
                ]
            )
        )

        if not in_tiered_vcf:
            continue

        info_out: Dict[str, object] = dict(c.info)
        info_out.update(
            {
                "TIER": tier_label,
                "TYPE": vtype,
                "DP": dp,
                "AD": f"{ad_ref},{ad_alt}",
                "VAF": f"{vaf:.5f}",
                "ALT_F": alt_f,
                "ALT_R": alt_r,
                "REF_F": int(counts["REF_F"]),
                "REF_R": int(counts["REF_R"]),
                "SB": f"{sb:.3g}" if sb is not None else ".",
                "CTX": ctx if ctx else ".",
                "APOBEC3": "Yes" if is_apobec else "No",
            }
        )

        vcf_body_lines.append(
            "\t".join(
                [
                    c.chrom,
                    str(c.pos),
                    ".",
                    c.ref,
                    c.alt,
                    qual_out,
                    "PASS",
                    _fmt_info(info_out),
                    fmt,
                    smp,
                ]
            )
        )

        rows_out.append(
            {
                "POS": c.pos,
                "REF": c.ref,
                "ALT": c.alt,
                "TYPE": vtype,
                "TIER": tier_label,
                "QUAL": qual_out,
                "DP": dp,
                "AD": f"{ad_ref},{ad_alt}",
                "VAF": f"{vaf:.5f}",
                "SB": f"{sb:.3g}" if sb is not None else "NA",
                "CTX": ctx or "",
                "APOBEC3": "Yes" if is_apobec else "No",
            }
        )

    def vcf_sort_key(line: str) -> Tuple[int, str, str]:
        p = line.split("\t")
        pos = int(p[1])
        ref = p[3]
        alt = p[4]
        return (pos, ref, alt)

    vcf_body_lines.sort(key=vcf_sort_key)
    unf_body_lines.sort(key=vcf_sort_key)
    rows_out.sort(key=lambda rr: (int(rr["POS"]), str(rr["REF"]), str(rr["ALT"])))

    out_unf = Path(f"{prefix}.variants.unfiltered.vcf")
    with open(out_unf, "wt", encoding="utf-8") as out:
        for hl in merged_header_unf:
            out.write(hl + "\n")
        for line in unf_body_lines:
            out.write(line + "\n")

    out_vcf = Path(f"{prefix}.variants.filtered.vcf")
    with open(out_vcf, "wt", encoding="utf-8") as out:
        for hl in merged_header:
            out.write(hl + "\n")
        for line in vcf_body_lines:
            out.write(line + "\n")

    out_html = Path(f"{prefix}.variants.html")

    settings = {
        "sample": args.sample,
        "virus": args.virus,
        "uniform_filters": {
            "min_qual": args.min_qual,
            "min_bq": args.min_bq,
            "min_mq": args.min_mq,
            "min_dp": args.min_dp,
            "min_alt_reads": args.min_alt_reads,
            "major_vaf": args.major_vaf,
            "minor_vaf_min": args.minor_vaf_min,
            "minor_vaf_max": args.minor_vaf_max,
            "min_sb_pvalue": args.min_sb_pvalue,
        },
        "inputs": {
            "bam": str(args.bam),
            "ref_fasta": str(args.ref_fasta),
            "clair3_vcf": str(args.clair3_vcf),
        },
        "reference": {"contig": contig, "length": contig_len},
    }
    _build_html_report(sample=args.sample, virus=args.virus, settings=settings, rows=rows_out, out_path=out_html)

    bam.close()
    fasta.close()


if __name__ == "__main__":
    main()

