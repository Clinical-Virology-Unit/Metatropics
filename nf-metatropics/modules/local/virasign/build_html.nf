process VIRASIGN_SUMMARY {
    tag "Metatropics summary report"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'library://jansendaan94_v2/metatropics/virasign:latest':
        'daanjansen94/virasign:latest' }"

    // Bind outdir and read results via params.
    def outAbs = file(params.outdir).toAbsolutePath().toString()

    // Ensure any control/blind files are visible in-container.
    def extraBindDirs = []
    def addBindDir = { String p ->
        if (!p) return
        def fp = file(p)
        def d = fp.isDirectory() ? fp : fp.parent
        if (d != null) extraBindDirs << d.toAbsolutePath().toString()
    }
    def zc = params.virasign_zscore_controls?.toString()?.trim()
    if (zc) {
        zc.split(',').collect { it.trim() }.findAll { it }.each { addBindDir(it) }
    }
    addBindDir(params.virasign_blind?.toString()?.trim())
    extraBindDirs = extraBindDirs.unique()

    if (workflow.containerEngine == 'docker' || workflow.containerEngine == 'podman') {
        def extra = extraBindDirs.collect { " -v ${it}:${it}" }.join('')
        containerOptions "-v ${outAbs}:${outAbs}${extra}"
    }
    if (workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer') {
        def extra = extraBindDirs.collect { " --bind ${it}:${it}" }.join('')
        containerOptions "--bind ${outAbs}:${outAbs}${extra}"
    }

    input:
    val _virasign_final_json_count

    output:
    path "Metatropics_Summary_*.html", emit: html
    path "Metatropics_Summary_*.csv",  emit: csv
    path "versions.yml",           emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def rawDbArg = params.virasign_database?.toString()?.trim()
    def effectiveDbArg = rawDbArg ?: 'RVDB'
    def virasignDbLabel = effectiveDbArg.replaceAll(/[^A-Za-z0-9._-]+/, '_')
    def outroot = file("${params.outdir}/Classification/virasign/${virasignDbLabel}").toAbsolutePath().toString()
    def outHtml = "Metatropics_Summary_${virasignDbLabel}.html"
    def outCsv  = "Metatropics_Summary_${virasignDbLabel}.csv"
    def zscore = params.virasign_zscore ?: 'true'
    def zscore_controls = params.virasign_zscore_controls?.toString()?.trim()
    def zscore_controls_safe = zscore_controls ?: ''
    def virasign_blind_safe = params.virasign_blind?.toString()?.trim() ?: ''
    def zscore_json = groovy.json.JsonOutput.toJson(zscore?.toString()?.trim() ?: 'false')
    def zscore_controls_json = groovy.json.JsonOutput.toJson(zscore_controls_safe)
    def virasign_blind_json = groovy.json.JsonOutput.toJson(virasign_blind_safe)
    def quality = params.quality?.toString() ?: 'NA'
    def depth = params.depth?.toString() ?: 'NA'
    def agreement = params.agreement?.toString() ?: 'NA'
    def opt = []
    def add = { c, f -> if (c) opt << f }
    add(!!zscore, "--zscore ${zscore}")
    add(zscore_controls, "--zscore-controls '${zscore_controls}'")
    def tail = task.ext.args?.toString()?.trim()
    def cmd = (['virasign', '--build-html', '-o', "${outroot}"] + opt + (tail ? [tail] : [])).join(' ')

    """
    OUT='${outroot}'
    # Short wait for slow FS after the last classification merge.
    SEC=0
    while [ "\$SEC" -lt 60 ]; do
      if [ -d "\$OUT" ] && find "\$OUT" -name '*_final_selected_references.json' \\( -type f -o -type l \\) 2>/dev/null | head -n1 | grep -q .; then
        break
      fi
      sleep 2
      SEC=\$((SEC + 2))
    done
    if [ ! -d "\$OUT" ]; then
      echo "ERROR: Virasign results root missing: \$OUT" >&2
      exit 1
    fi
    if ! find "\$OUT" -name '*_final_selected_references.json' \\( -type f -o -type l \\) 2>/dev/null | head -n1 | grep -q .; then
      # No confident hits in any sample: don't fail the pipeline; emit a small placeholder summary.
      csv="${outCsv}"
      html="${outHtml}"
      printf "sample,confident\\n" > "\$csv"
      cat > "\$html" <<'EOF'
<!doctype html>
<html lang="en">
  <head>
    <meta charset="utf-8"/>
    <meta name="viewport" content="width=device-width, initial-scale=1"/>
    <title>Virasign summary</title>
  </head>
  <body>
    <h1>Virasign summary</h1>
    <p>No confident hits were found (no <code>*_final_selected_references.json</code> files present).</p>
  </body>
</html>
EOF
      cat > versions.yml <<'END_VERSIONS'
"${task.process}":
    virasign: \$(virasign --version 2>&1 | head -n1 || echo 'unknown')
END_VERSIONS
      exit 0
    fi

    ${cmd}
    cp -f "${outroot}"/results_summary_*.html ./__results_summary.html
    cp -f "${outroot}"/results_summary_*.csv  ./__results_summary.csv
    mv -f ./__results_summary.html "${outHtml}"
    mv -f ./__results_summary.csv  "${outCsv}"

# Compute consensus-derived breadth from draft consensus FASTA and merge into summary CSV.
    python3 - <<'PY'
import csv
import json
import gzip
import os
import re
from pathlib import Path

outdir = Path("${params.outdir}")
summary_csv = Path("${outCsv}")
summary_html = Path("${outHtml}")
readcount_csv = outdir / "Summary" / "readcount" / "read_counts.csv"
cons_root = outdir / "Consensus" / "bcftools"

virasign_runtime_log = outdir / "Classification" / "virasign" / "${virasignDbLabel}" / ".virasign.log"

# Virasign background/z-score configuration (for human-readable documentation in HTML).
zscore_enabled_raw = ${zscore_json}
zscore_controls_text = ${zscore_controls_json}
virasign_blind_text = ${virasign_blind_json}
zscore_enabled = str(zscore_enabled_raw).strip().lower() in {"true","1","yes","y"}
zscore_controls_text = str(zscore_controls_text or "").strip()
virasign_blind_text = str(virasign_blind_text or "").strip()

if not zscore_enabled:
    zscore_water_msg = f"Z-score disabled (--zscore={zscore_enabled_raw})."
else:
    if zscore_controls_text:
        zscore_water_msg = f"Water controls (explicit --zscore-controls): {zscore_controls_text}"
    else:
        zscore_water_msg = "Water controls: auto-detected by Virasign (no --zscore-controls provided)."

if virasign_blind_text:
    blind_msg = f"Blind/background file (explicit --virasign_blind): {virasign_blind_text}"
else:
    blind_msg = "Blind/background file: not provided (Virasign uses defaults/auto when needed)."

# If Virasign couldn't compute Z-scores due to missing/insufficient water controls,
# the produced summary CSV will have an empty Z-score column.
# Detect that so the HTML documentation can explicitly say "none detected/used".
zscore_has_values = False
try:
    if summary_csv.exists():
        with summary_csv.open("r", newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            if reader.fieldnames and ("Z-score" in reader.fieldnames or "Z-score " in reader.fieldnames):
                zscore_col = "Z-score" if "Z-score" in reader.fieldnames else "Z-score "
                for row in reader:
                    raw_v = (row.get(zscore_col) or "").strip()
                    if raw_v and raw_v not in {"—", "-"}:
                        zscore_has_values = True
                        break
except Exception:
    pass

if zscore_enabled and not zscore_has_values:
    if zscore_controls_text:
        zscore_water_msg = (
            f"Water controls provided, but none were usable/available for Z-score computation "
            f"(Virasign left Z-score empty). Provided: {zscore_controls_text}"
        )
    else:
        zscore_water_msg = (
            "Water controls: none usable/available (Virasign skipped Z-score computation; Z-score column is empty)."
        )

# Values for CSV annotation (computed here after inspecting outputs)
background_used = "yes" if (zscore_enabled and zscore_has_values) else "no"

# Determine which samples are background controls for Z-score.
# Prefer explicit `--zscore-controls` (paths or sample IDs). If not provided,
# fall back to common naming patterns (H2O/water/blank/control).
def _norm_sample_id(s: str) -> str:
    s = (s or "").strip()
    if not s:
        return ""
    for suf in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
        if s.endswith(suf):
            s = s[: -len(suf)]
            break
    # NOTE: Use \\Z instead of an end-of-line anchor to avoid Groovy interpolation issues.
    s = re.sub(r"\\.fastp\\Z", "", s)
    s = re.sub(r"(_T\\d+_other_T\\d+|_other_T\\d+|_T\\d+_other|_T\\d+|_other)\\Z", "", s)
    return s

def _parse_controls_list(txt: str) -> list:
    if not txt:
        return []
    out = []
    for part in str(txt).split(","):
        p = part.strip()
        if not p:
            continue
        p = os.path.basename(p)
        out.append(_norm_sample_id(p))
    return [x for x in out if x]

controls_norm = set(_parse_controls_list(zscore_controls_text))

def is_background_sample(sample: str) -> bool:
    s_norm = _norm_sample_id(sample)
    if not s_norm:
        return False
    if controls_norm:
        return s_norm in controls_norm
    s = s_norm.lower()
    # Heuristic fallback: be strict and only match typical water/blank controls.
    # Do NOT match generic substrings like "neg" because many studies encode that in sample names
    # (e.g. LASVnegRUN1_...) and it would incorrectly flag all samples as background.
    return any(k in s for k in ("h2o", "water", "blank"))

# For auto-detection mode, show exact Virasign log lines where available.
# This can include explicit control files selected or "none/insufficient controls".
zscore_runtime_lines = []
try:
    if virasign_runtime_log.exists():
        seen = set()
        with virasign_runtime_log.open("r", encoding="utf-8", errors="ignore") as handle:
            for line in handle:
                s = line.strip()
                if not s:
                    continue
                lo = s.lower()
                if ("z-score" in lo) or ("zscore" in lo) or ("water control" in lo) or ("zscore-controls" in lo):
                    if s not in seen:
                        seen.add(s)
                        zscore_runtime_lines.append(s)
except Exception:
    pass

if zscore_enabled and not zscore_controls_text and zscore_runtime_lines:
    zscore_water_msg = "Water controls: auto-detected by Virasign. " + " | ".join(zscore_runtime_lines[:3])

def safe_pct(num, den):
    if den == 0:
        return ""
    return f"{(100.0 * num / den):.2f}"

def parse_fasta_counts(path: Path):
    total = 0
    called = 0
    n_bases = 0
    gap_bases = 0
    acgt = 0
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line or line.startswith(">"):
                continue
            seq = line.strip().upper()
            for ch in seq:
                if ch.isspace():
                    continue
                total += 1
                if ch == "N":
                    n_bases += 1
                else:
                    called += 1
                if ch == "-":
                    gap_bases += 1
                if ch in {"A", "C", "G", "T"}:
                    acgt += 1
    return total, called, n_bases, gap_bases, acgt

vsign_parent = outdir / "Classification" / "virasign"

def sanitize_species_slug(raw):
    if raw is None:
        return ""
    s = str(raw).strip()
    if not s:
        return ""
    s = re.sub(r"[^A-Za-z0-9._-]+", "_", s)
    s = re.sub(r"_+", "_", s).strip("_")
    return s

def virus_slug_from_hit(hit):
    acc = str(hit.get("accession") or "").strip()
    if not acc:
        return ""
    raw_sp = hit.get("organism") or hit.get("viral_species") or ""
    raw_sp = str(raw_sp).strip() if raw_sp else ""
    if not raw_sp and hit.get("description"):
        raw_sp = str(hit.get("description")).strip()[:120]
    sp = sanitize_species_slug(raw_sp)
    return f"{acc}_{sp}" if sp else acc

def load_slug_to_accession():
    m = {}
    if not vsign_parent.is_dir():
        return m
    for js in sorted(vsign_parent.rglob("*_final_selected_references.json")):
        try:
            data = json.loads(js.read_text(encoding="utf-8"))
        except Exception:
            continue
        if not isinstance(data, list):
            continue
        for hit in data:
            if not isinstance(hit, dict):
                continue
            acc = str(hit.get("accession") or "").strip()
            if not acc:
                continue
            slug = virus_slug_from_hit(hit)
            if slug:
                m[slug] = acc
    return m

slug_to_accession = load_slug_to_accession()

rows = []
by_pair = {}
if cons_root.exists():
    for fasta in sorted(cons_root.rglob("*.consensus.fasta")):
        sample = fasta.parent.name
        base = fasta.name
        suffix = ".consensus.fasta"
        virus_key = ""
        if base.endswith(suffix):
            core = base[:-len(suffix)]
            prefix = f"{sample}."
            if core.startswith(prefix):
                virus_key = core[len(prefix):]
            else:
                parts = core.split(".", 1)
                virus_key = parts[1] if len(parts) == 2 else core

        bare_acc = slug_to_accession.get(virus_key, virus_key)
        total, called, n_bases, gap_bases, acgt = parse_fasta_counts(fasta)

        # Prefer reference-coordinate breadth computed during consensus building.
        # This avoids >100% due to insertions and counts deletions as determined.
        breadth = ""
        metrics = fasta.with_suffix(".metrics.json")
        if metrics.exists():
            try:
                m = json.loads(metrics.read_text(encoding="utf-8"))
                v = m.get("consensus_breadth_pct")
                if v is not None and str(v).strip() != "":
                    breadth = f"{float(v):.2f}"
            except Exception:
                breadth = ""

        # Fallback: estimate unambiguous breadth from the FASTA itself (A/C/G/T only).
        # Note: this can be impacted by indels, so it is used only when metrics are missing.
        if breadth == "":
            denom = total - gap_bases
            breadth = safe_pct(acgt, denom)
        strict_breadth = breadth
        rec = {
            "sample": sample,
            "accession": bare_acc,
            "consensus_breadth_pct": breadth,
            "consensus_acgt_breadth_pct": strict_breadth,
        }
        rows.append(rec)
        by_pair[(sample, bare_acc)] = rec
        if virus_key and virus_key != bare_acc:
            by_pair[(sample, virus_key)] = rec

qc_reads_by_sample = {}
readlen_med_by_sample = {}

# For the Metatropics Virasign report, "QC reads" should always mean:
# the number of reads that were actually used as INPUT for Virasign
# (post-fastplong, and post any enabled host depletion).
def strip_extensions(name: str) -> str:
    while True:
        # NOTE: Use \\Z end-anchor to avoid Groovy interpolation inside Nextflow strings.
        new = re.sub(r"\\.(fastq|fq|fastp|gz|meta|csv)\\Z", "", name)
        if new == name:
            return name
        name = new

def extract_sample_name(filename: str) -> str:
    name = strip_extensions(filename)
    name = re.sub(r"_(human_depleted|host_depleted)\\Z", "", name)
    name = re.sub(r"_(human|host)\\Z", "", name)
    name = re.sub(r"_other\\Z", "", name)
    name = re.sub(r"_viral\\Z", "", name)
    name = re.sub(r"_classification_results\\Z", "", name)
    name = re.sub(r"_fixed\\Z", "", name)
    name = re.sub(r"\\.fastp\\Z", "", name)
    return name

def count_fastq_gz_reads(path: Path) -> int:
    n = 0
    with gzip.open(path, "rt", errors="replace") as fh:
        while True:
            l1 = fh.readline()
            if not l1:
                break
            fh.readline()
            fh.readline()
            fh.readline()
            n += 1
    return n

def median_from_counts(counts: dict) -> str:
    # Exact median from a {length: count} histogram.
    if not counts:
        return ""
    total = sum(int(v) for v in counts.values() if v)
    if total <= 0:
        return ""
    # 0-based median indices in the sorted multiset
    lo = (total - 1) // 2
    hi = total // 2
    seen = 0
    lo_val = None
    hi_val = None
    for length in sorted(counts.keys()):
        c = int(counts.get(length) or 0)
        if c <= 0:
            continue
        nxt = seen + c
        if lo_val is None and lo < nxt:
            lo_val = length
        if hi_val is None and hi < nxt:
            hi_val = length
        if lo_val is not None and hi_val is not None:
            break
        seen = nxt
    if lo_val is None or hi_val is None:
        return ""
    return str(int(round((lo_val + hi_val) / 2.0)))

def fastq_gz_length_counts(path: Path) -> dict:
    # Exact scan over all reads; store only a histogram (length -> count).
    counts = {}
    with gzip.open(path, "rt", errors="replace") as fh:
        while True:
            h = fh.readline()
            if not h:
                break
            seq = fh.readline()
            if not seq:
                break
            fh.readline()
            fh.readline()
            L = len(seq.strip())
            counts[L] = counts.get(L, 0) + 1
    return counts

def count_dir(directory: Path, pattern: re.Pattern) -> dict:
    out = {}
    length_counts = {}
    if not directory.exists():
        return out
    for p in directory.iterdir():
        if not p.is_file():
            continue
        if not pattern.search(p.name):
            continue
        sample = extract_sample_name(p.name)
        out[sample] = out.get(sample, 0) + count_fastq_gz_reads(p)
        try:
            c = fastq_gz_length_counts(p)
            if c:
                agg = length_counts.setdefault(sample, {})
                for k, v in c.items():
                    agg[k] = agg.get(k, 0) + int(v)
        except Exception:
            pass
    # Store per-sample median read length.
    for s, counts in length_counts.items():
        readlen_med_by_sample[s] = median_from_counts(counts)
    return out

reads_root = outdir / "Reads"
# Prefer the most-depleted stage if present (this matches readsForViralClassification).
qc_reads_by_sample = count_dir(reads_root / "nohost", re.compile(r"_(host_depleted|other)\\.fastq\\.gz\\Z"))
if not qc_reads_by_sample:
    qc_reads_by_sample = count_dir(reads_root / "nohuman", re.compile(r"_(human_depleted|other)\\.fastq\\.gz\\Z"))
if not qc_reads_by_sample:
    qc_reads_by_sample = count_dir(reads_root / "fastplong", re.compile(r"\\.fastp\\.fastq\\.gz\\Z"))

# As a last resort, fall back to whatever readcount produced (if available).
if not qc_reads_by_sample and readcount_csv.exists():
    try:
        with readcount_csv.open("r", newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            for row in reader:
                sample = (row.get("sample") or "").strip()
                v = (row.get("trimmed_reads") or "").strip()
                if sample:
                    qc_reads_by_sample[sample] = v
    except Exception:
        qc_reads_by_sample = {}

qc_map_json = json.dumps(qc_reads_by_sample, separators=(",", ":"))
readlen_map_json = json.dumps(readlen_med_by_sample, separators=(",", ":"))

if summary_csv.exists():
    with summary_csv.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        in_fields = list(reader.fieldnames or [])
        data = list(reader)

    low = {f.lower().strip(): f for f in in_fields}
    sample_col = None
    for key in ("sample", "sample_id", "sampleid"):
        if key in low:
            sample_col = low[key]
            break
    accession_col = None
    for key in ("accession", "virus", "reference", "reference_id", "ref"):
        if key in low:
            accession_col = low[key]
            break

    qc_header = "QC reads"
    readlen_header = "Read Length (Med)"
    consensus_header = "Consensus Breadth (%)"
    bg_used_header = "Background"

    out_fields = list(in_fields)

    # Remove any old technical consensus header so we can re-add consistently.
    out_fields = [f for f in out_fields if f != "consensus_breadth_pct"]
    out_fields = [f for f in out_fields if f != bg_used_header]

    # Place Background after Sample when Sample column exists, otherwise keep first.
    if sample_col and sample_col in out_fields:
        sample_idx = out_fields.index(sample_col)
        out_fields.insert(sample_idx + 1, bg_used_header)
    else:
        out_fields = [bg_used_header] + out_fields

    # Enforce QC column before mapped reads.
    if qc_header in out_fields:
        out_fields.remove(qc_header)
    if "Mapped Reads (#)" in out_fields:
        mapped_idx = out_fields.index("Mapped Reads (#)")
        out_fields.insert(mapped_idx, qc_header)
    else:
        out_fields.append(qc_header)

    # Insert median read length after Avg Identity (%) when present.
    if readlen_header in out_fields:
        out_fields.remove(readlen_header)
    if "Avg Identity (%)" in out_fields:
        idx = out_fields.index("Avg Identity (%)")
        out_fields.insert(idx + 1, readlen_header)
    else:
        out_fields.append(readlen_header)

    # Add consensus column as the last column.
    if consensus_header in out_fields:
        out_fields.remove(consensus_header)
    out_fields.append(consensus_header)

    if sample_col:
        # Optional sample-level fallback when each sample has a single consensus.
        by_sample = {}
        for rec in rows:
            by_sample.setdefault(rec["sample"], []).append(rec)
        for row in data:
            sample = (row.get(sample_col) or "").strip()
            acc = (row.get(accession_col) or "").strip() if accession_col else ""
            hit = by_pair.get((sample, acc))
            if hit is None and not acc:
                sample_hits = by_sample.get(sample, [])
                if len(sample_hits) == 1:
                    hit = sample_hits[0]
            if zscore_enabled and zscore_has_values:
                row[bg_used_header] = "yes" if is_background_sample(sample) else "no"
            else:
                # If Z-score is not used (disabled or unusable), make this explicit instead of leaving blank.
                row[bg_used_header] = "no"
            row[qc_header] = qc_reads_by_sample.get(sample, "")
            row[readlen_header] = readlen_med_by_sample.get(sample, "")
            row[consensus_header] = "" if hit is None else hit["consensus_breadth_pct"]
    else:
        for row in data:
            row[qc_header] = ""
            row[bg_used_header] = "no"
            row[readlen_header] = ""
            row[consensus_header] = ""

    with summary_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=out_fields)
        writer.writeheader()
        writer.writerows(data)

if summary_html.exists():
    import json

    text = summary_html.read_text(encoding="utf-8")

    # If the pipeline is re-run, remove any previously injected consensus panel.
    # Use plain string operations (avoid regex backslash parsing issues).
    section_start = text.find('<section id="consensus-breadth"')
    if section_start != -1:
        section_end = text.find('</section>', section_start)
        if section_end != -1:
            section_end = section_end + len('</section>')
            # Often the section is immediately followed by a <script>...</script>.
            script_start = text.find('<script', section_end)
            if script_start != -1:
                script_end = text.find('</script>', script_start)
                if script_end != -1:
                    script_end = script_end + len('</script>')
                    text = text[:section_start] + text[script_end:]
                else:
                    text = text[:section_start] + text[section_end:]
            else:
                text = text[:section_start] + text[section_end:]

    # Ensure the per-sample "Download Table (CSV)" header matches the appended column.
    text = text.replace(
        "let csv = 'Accession,Organism,Viral Species,Nextclade Clade,Segment,Mapped Reads (#),Identity (%),Coverage Depth (x),Coverage Breadth (%),NOGR (#/bases),Z-score\\n';",
        "let csv = 'Accession,Organism,Viral Species,Nextclade Clade,Segment,Mapped Reads (#),Identity (%),Coverage Depth (x),Coverage Breadth (%),NOGR (#/bases),Z-score,Consensus Breadth (%)\\n';",
        1
    )
    if "Consensus Breadth (%)" not in text:
        text = text.replace("Z-score\\n';", "Z-score,Consensus Breadth (%)\\n';", 1)

    # Compact map used by client-side JS to populate the new table column.
    # Structure: { "<sample>": { "<accession>": <float breadth_pct> } }
    cons_map = {}
    for rec in rows:
        s = (rec.get("sample") or "").strip()
        a = (rec.get("accession") or "").strip()
        v = rec.get("consensus_breadth_pct")
        if not s or not a or v in (None, ""):
            continue
        try:
            cons_map.setdefault(s, {})[a] = float(v)
        except Exception:
            pass
    cons_map_json = json.dumps(cons_map, separators=(",", ":"))

    # Sample -> Background (yes/no) map for CSV downloads triggered from HTML.
    bg_map = {}
    try:
        if sample_col:
            for row in data:
                s = (row.get(sample_col) or "").strip()
                if not s:
                    continue
                bg_map[s] = (row.get(bg_used_header) or "").strip()
    except Exception:
        bg_map = {}
    bg_map_json = json.dumps(bg_map, separators=(",", ":"))

    inj_script = '''
<script>
(function(){
  const consensusMap = __CONS_MAP_JSON__;
  const qcReadsMap = __QC_READS_MAP_JSON__;
  const readLenMap = __READLEN_MAP_JSON__;
  const bgMap = __BG_MAP_JSON__;

  function getVal(sampleName, accession){
    if (!consensusMap || !consensusMap[sampleName]) return null;
    const v = consensusMap[sampleName][accession];
    if (v === undefined || v === null || v === '') return null;
    return Number(v);
  }

  function getQcReads(sampleName){
    if (!qcReadsMap) return '';
    const v = qcReadsMap[sampleName];
    return (v === undefined || v === null) ? '' : String(v);
  }

  function getReadLenMed(sampleName){
    if (!readLenMap) return '';
    const v = readLenMap[sampleName];
    return (v === undefined || v === null) ? '' : String(v);
  }

  function getBackground(sampleName){
    if (!bgMap) return '';
    const v = bgMap[sampleName];
    return (v === undefined || v === null) ? '' : String(v);
  }

  function ensureHeader(table){
    if (!table) return;
    const headerRows = table.querySelectorAll('thead tr');
    if (headerRows.length < 2) return;
    const filterRow = headerRows[0];
    const mainRow = headerRows[1];
    const colgroup = table.querySelector('colgroup');

    if (colgroup && !colgroup.querySelector('col[data-col="consensus_breadth_pct"]')) {
      const col = document.createElement('col');
      col.setAttribute('data-col', 'consensus_breadth_pct');
      col.style.width = '140px';
      col.style.display = 'none';
      colgroup.appendChild(col);
    }

    if (filterRow && !filterRow.querySelector('th[data-col="consensus_breadth_pct"]')) {
      const filterTh = document.createElement('th');
      filterTh.setAttribute('data-col', 'consensus_breadth_pct');
      filterTh.style.display = 'none';

      // Add a numeric filter for consensus breadth when the column is toggled on.
      const inp = document.createElement('input');
      inp.type = 'text';
      inp.placeholder = 'Min Breadth';
      inp.addEventListener('keyup', function () {
        // applyAllFilters is defined in the base HTML (already on window)
        if (typeof window.applyAllFilters === 'function') window.applyAllFilters(sampleName);
      });
      inp.style.width = '100%';
      inp.style.padding = '5px 6px';
      inp.style.border = '1px solid #ddd';
      inp.style.borderRadius = '4px';
      inp.style.fontSize = '0.9em';
      inp.style.boxSizing = 'border-box';
      filterTh.appendChild(inp);

      filterRow.appendChild(filterTh);
    }

    if (mainRow.querySelector('th[data-col="consensus_breadth_pct"]')) return;

    const th = document.createElement('th');
    th.textContent = 'Consensus Breadth (%)';
    th.setAttribute('data-col', 'consensus_breadth_pct');
    th.className = 'stats sortable';
    th.setAttribute('data-sort', 'consensus_breadth_pct');
    th.style.display = 'none';
    mainRow.appendChild(th);
  }

  function ensureToggleButton(sampleName){
    const sampleSection = document.getElementById('sample-' + sampleName);
    if (!sampleSection) return;
    const header = sampleSection.querySelector('.table-header');
    if (!header) return;
    if (header.querySelector('button[data-role="consensus-metrics-toggle"]')) return;

    const button = document.createElement('button');
    button.type = 'button';
    button.setAttribute('data-role', 'consensus-metrics-toggle');
    button.textContent = 'Consensus Metrics';
    button.addEventListener('click', function(){
      toggleConsensusColumn(sampleName);
    });

    // Keep download button + consensus toggle grouped on the right.
    let rightActions = header.querySelector('div[data-role="consensus-right-actions"]');
    const downloadBtn = header.querySelector('button[onclick^="downloadTableAsCSV"]');
    if (!rightActions) {
      rightActions = document.createElement('div');
      rightActions.setAttribute('data-role', 'consensus-right-actions');
      rightActions.style.display = 'flex';
      rightActions.style.gap = '0.6rem';
      rightActions.style.alignItems = 'center';

      if (downloadBtn) {
        header.insertBefore(rightActions, downloadBtn);
        rightActions.appendChild(downloadBtn);
      } else {
        header.appendChild(rightActions);
      }
    }

    rightActions.appendChild(button);
  }

  function setConsensusColumnVisibility(sampleName, visible){
    const table = document.getElementById('table-' + sampleName);
    const tbody = document.getElementById('tbody-' + sampleName);
    if (!table) return;

    const headerRows = table.querySelectorAll('thead tr');
    if (headerRows.length < 2) return;
    const filterRow = headerRows[0];
    const mainRow = headerRows[1];

    // 1) Consensus column itself
    table.querySelectorAll('th[data-col="consensus_breadth_pct"], td[data-col="consensus_breadth_pct"]').forEach(el => {
      el.style.display = visible ? '' : 'none';
    });

    // 1b) Spacing fix: when table-layout is fixed, hiding only th/td can leave
    // whitespace reserved for hidden columns. Hide the <col> elements too.
    const colgroup = table.querySelector('colgroup');
    if (colgroup) {
      const cols = Array.from(colgroup.querySelectorAll('col'));
      if (!visible) {
        // Default view: keep original layout, but do NOT reserve width for the
        // consensus column (we already hide the corresponding th/td above).
        table.style.tableLayout = 'fixed';
        cols.forEach(c => {
          if (c.getAttribute('data-col') === 'consensus_breadth_pct') {
            c.style.display = 'none';
          } else {
            c.style.display = '';
          }
        });
      } else {
        table.style.tableLayout = 'auto';

        const ths = Array.from(mainRow.querySelectorAll('th'));
        const n = Math.min(cols.length, ths.length);
        for (let i = 0; i < n; i++) {
          const th = ths[i];
          if (!th) {
            cols[i].style.display = 'none';
            continue;
          }
          const dataCol = th.getAttribute('data-col');
          const dataSort = th.getAttribute('data-sort');
          const show =
            (dataCol === 'consensus_breadth_pct') ||
            (dataSort === 'breadth') ||
            (dataSort === 'nogr_regions') ||
            (dataSort === 'zscore') ||
            (!dataSort); // base columns have no data-sort
          cols[i].style.display = show ? '' : 'none';
        }
        // If the colgroup has extra columns, hide them.
        for (let i = n; i < cols.length; i++) {
          cols[i].style.display = 'none';
        }
      }
    }

    // 2) Hide/show the other columns when consensus view is enabled
    //    Compact view (keep these):
    //      Accession | Organism | Viral species | Clade | Segment | Breadth | Consensus Breadth | NOGR | Z-score
    //    Hide only: Mapped reads (#), Identity, Depth.
    const hideMainSorts = ['mapped_reads','identity','depth'];
    hideMainSorts.forEach(s => {
      const th = mainRow.querySelector('th[data-sort="' + s + '"]');
      if (th) th.style.display = visible ? 'none' : '';
    });

    const hidePlaceholders = new Set(['Min Reads','Min ID','Min Depth']);
    if (filterRow) {
      filterRow.querySelectorAll('input[type="text"]').forEach(inp => {
        const ph = inp.getAttribute('placeholder') || '';
        const th = inp.closest('th');
        if (!th) return;
        if (hidePlaceholders.has(ph)) th.style.display = visible ? 'none' : '';
      });
    }

    if (tbody) {
      tbody.querySelectorAll('tr').forEach(row => {
        const tds = row.querySelectorAll('td');
        // Expected base columns (0-based):
        // 0 Accession, 1 Organism, 2 Viral Species, 3 Clade, 4 Segment,
        // 5 Mapped Reads, 6 Identity, 7 Depth, 8 Breadth,
        // 9 NOGR, 10 Z-score, 11 Consensus Breadth (appended by this script)
        if (tds.length > 10) {
          // Hide columns for compact view (visible=true)
          tds[5].style.display = visible ? 'none' : '';
          tds[6].style.display = visible ? 'none' : '';
          tds[7].style.display = visible ? 'none' : '';
          // Keep NOGR and Z-score visible
          tds[9].style.display = '';
          tds[10].style.display = '';
        }
      });
    }

    const sampleSection = document.getElementById('sample-' + sampleName);
    if (sampleSection) {
      const button = sampleSection.querySelector('button[data-role="consensus-metrics-toggle"]');
      if (button) {
        button.textContent = visible ? 'Hide Consensus Metrics' : 'Consensus Metrics';
      }
    }
  }

  function toggleConsensusColumn(sampleName){
    const table = document.getElementById('table-' + sampleName);
    if (!table) return;
    const header = table.querySelector('th[data-col="consensus_breadth_pct"]');
    const visible = !!(header && header.style.display !== 'none');
    setConsensusColumnVisibility(sampleName, !visible);
  }

  function injectForSample(sampleName){
    const table = document.getElementById('table-' + sampleName);
    const tbody = document.getElementById('tbody-' + sampleName);
    if (!table || !tbody) return;

    ensureHeader(table);
    ensureToggleButton(sampleName);

    tbody.querySelectorAll('tr').forEach(row => {
      const accession = row.getAttribute('data-accession') || '';
      const v = getVal(sampleName, accession);
      const txt = (v === null) ? '-' : (v.toFixed(2) + '%');

      let td = row.querySelector('td[data-col="consensus_breadth_pct"]');
      if (!td){
        td = document.createElement('td');
        td.setAttribute('data-col','consensus_breadth_pct');
        td.className = 'stats';
        td.style.display = 'none';
        row.appendChild(td);
      }
      row.setAttribute('data-consensus-breadth', v === null ? '' : String(v));
      td.textContent = txt;
    });

    setConsensusColumnVisibility(sampleName, false);
  }

  const origPopulateTable = window.populateTable;
  if (typeof origPopulateTable === 'function'){
    window.populateTable = function(sampleName){
      origPopulateTable(sampleName);
      injectForSample(sampleName);
    };
  }

  // (No Z-score background metadata injection: keep HTML clean.)

  // Override downloadTableAsCSV so the downloaded CSV includes QC reads + median read length + consensus breadth.
  const origDownloadTableAsCSV = window.downloadTableAsCSV;
  window.downloadTableAsCSV = function(sampleName) {
    const tbody = document.getElementById('tbody-' + sampleName);
    if (!tbody) return;

    const NL = String.fromCharCode(10);
    let csv = 'Sample,Background,Accession,Organism,Viral Species,Nextclade Clade,Segment,QC reads,Mapped Reads (#),Identity (%),Read Length (Med),Coverage Depth (x),Coverage Breadth (%),NOGR (#/bases),Z-score,Consensus Breadth (%)' + NL;
    const rows = tbody.querySelectorAll('tr:not([style*="display: none"])');

    rows.forEach(row => {
      const cells = row.querySelectorAll('td');
      const accession = row.getAttribute('data-accession') || '';
      const qc = getQcReads(sampleName);
      const rl = getReadLenMed(sampleName);
      const bg = getBackground(sampleName);
      const cons = getVal(sampleName, accession);
      const consTxt = (cons === null) ? '' : cons.toFixed(2);

      // Base table has 11 cells (up to Z-score). The injected consensus cell is appended at end
      // but may be hidden; we always compute it from consensusMap to keep download consistent.
      const out = [];
      for (let i = 0; i < Math.min(cells.length, 11); i++) {
        let text = cells[i].textContent.trim();
        // strip the chart icon used in breadth column
        text = text.replace(/📊/g, '').trim();
        if (text.includes(',') || text.includes('\"')) {
          text = '\"' + text.replace(/\"/g, '\"\"') + '\"';
        }
        out.push(text);
      }

      // Prepend Sample + Background, insert QC reads after Segment, append Consensus Breadth.
      out.unshift(sampleName);
      out.splice(1, 0, bg);
      out.splice(7, 0, qc);
      // After QC insertion, Identity (%) is at index 9 → insert read length at index 10.
      out.splice(10, 0, rl);
      out.push(consTxt);
      csv += out.join(',') + NL;
    });

    const blob = new Blob(['\\ufeff' + csv], { type: 'text/csv;charset=utf-8;' });
    const link = document.createElement('a');
    link.href = URL.createObjectURL(blob);
    link.download = sampleName + '_results.csv';
    link.click();
  };

  // Override the all-samples CSV export to include Background + QC reads + median read length + consensus breadth.
  // This does not change what is shown in the HTML tables; it only affects the downloaded file.
  if (typeof window.downloadAllSamplesCSV === 'function') {
    const origDownloadAll = window.downloadAllSamplesCSV;
    window.downloadAllSamplesCSV = function() {
      const NL = String.fromCharCode(10);
      let rows = [['Sample','Background','Accession','Organism','Viral Species','Nextclade Clade','Segment','QC reads','Mapped Reads (#)','Avg Identity (%)','Read Length (Med)','Coverage Depth (x)','Coverage Breadth (%)','NOGR (#/bases)','Z-score','Consensus Breadth (%)']];
      if (typeof allSamples === 'undefined' || typeof samplesData === 'undefined') return origDownloadAll();

      allSamples.forEach(sname => {
        const sdata = samplesData[sname];
        const qc = getQcReads(sname);
        const rl = getReadLenMed(sname);
        const bg = getBackground(sname);
        if (!sdata) { rows.push([sname,bg,'','','','','',qc,'','',rl,'','','','','']); return; }
        let refs = [];
        Object.keys(sdata).forEach(db => { if (sdata[db].references) refs = refs.concat(sdata[db].references); });
        if (refs.length === 0) { rows.push([sname,bg,'','','','','',qc,'','',rl,'','','','','']); return; }
        refs.sort((a,b) => (b.coverage_breadth||0)-(a.coverage_breadth||0));
        refs.forEach(r => {
          const acc = r.accession||'';
          const cons = getVal(sname, acc);
          rows.push([
            sname,
            bg,
            acc,
            r.organism||'',
            r.viral_species||'',
            r.nextclade_clade||'',
            r.segment||'',
            qc,
            r.mapped_reads||0,
            (r.avg_identity||0).toFixed(2),
            rl,
            (r.coverage_depth||0).toFixed(2),
            ((r.coverage_breadth||0)*100).toFixed(1),
            String((r.nogr_regions||r.non_overlapping_reads||0)) + '|' + String((r.nogr_bases||r.non_overlapping_bases||0)),
            (r.zscore !== undefined && r.zscore !== null) ? Number(r.zscore).toFixed(2) : '',
            (cons === null) ? '' : cons.toFixed(2)
          ]);
        });
      });
      const csv = rows.map(r => r.map(c => '\"'+String(c).replace(/\"/g,'\"\"')+'\"').join(',')).join(NL);
      const blob = new Blob([csv], {type:'text/csv'});
      const url = URL.createObjectURL(blob);
      const a = document.createElement('a');
      a.href = url;
      a.download = 'results_summary_all_samples.csv';
      a.click();
      URL.revokeObjectURL(url);
    };
  }

  // Extend filtering with an additional consensus-breadth rule, but only when
  // the consensus column + its filter input are visible.
  const origApplyAllFilters = window.applyAllFilters;
  if (typeof origApplyAllFilters === 'function'){
    window.applyAllFilters = function(sampleName){
      origApplyAllFilters(sampleName);
      try {
        const table = document.getElementById('table-' + sampleName);
        if (!table) return;
        const th = table.querySelector('th[data-col="consensus_breadth_pct"]');
        if (!th || th.style.display === 'none') return;
        const input = th.querySelector('input[type="text"]');
        if (!input) return;
        const raw = (input.value || '').trim();
        if (raw === '') return;
        const minVal = parseFloat(raw);
        if (isNaN(minVal)) return;

        const tbody = document.getElementById('tbody-' + sampleName);
        if (!tbody) return;

        tbody.querySelectorAll('tr').forEach(row => {
          // If some other filter already hid the row, don't re-show it here.
          if (row.style.display === 'none') return;
          const vRaw = row.getAttribute('data-consensus-breadth') || '';
          const v = parseFloat(vRaw);
          if (isNaN(v) || v < minVal) row.style.display = 'none';
        });
      } catch (e) {}
    };
  }

  // Inject into the currently active sample (if already rendered).
  try {
    const active = document.querySelector('.sample-section.active');
    if (active && active.id && active.id.startsWith('sample-')){
      const sn = active.id.slice('sample-'.length);
      injectForSample(sn);
    }
  } catch (e) {}
})();
</script>
'''

    inj_script = inj_script.replace("__CONS_MAP_JSON__", cons_map_json)
    inj_script = inj_script.replace("__QC_READS_MAP_JSON__", qc_map_json)
    inj_script = inj_script.replace("__READLEN_MAP_JSON__", readlen_map_json)
    inj_script = inj_script.replace("__BG_MAP_JSON__", bg_map_json)

    if 'data-col="consensus_breadth_pct"' not in text:
        if "</body>" in text:
            text = text.replace("</body>", inj_script + "\\n</body>", 1)
        else:
            text = text + "\\n" + inj_script

    summary_html.write_text(text, encoding="utf-8")
PY

    # Keep the shared Classification/virasign tree focused on per-sample outputs only.
    rm -f "${outroot}"/results_summary_*.html "${outroot}"/results_summary_*.csv "${outroot}"/.virasign.log || true

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        virasign: \$(virasign --version 2>&1 | head -n1 || echo 'unknown')
    END_VERSIONS
    """
}

