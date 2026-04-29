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

    # Compute consensus-derived breadth from polished FASTA and merge into summary CSV.
    python3 - <<'PY'
import csv
import os
from pathlib import Path

outdir = Path("${params.outdir}")
summary_csv = Path("${outCsv}")
summary_html = Path("${outHtml}")
readcount_csv = outdir / "Summary" / "readcount" / "read_counts.csv"
cons_root = outdir / "Consensus" / "homopolish"

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

rows = []
by_pair = {}
if cons_root.exists():
    for fasta in sorted(cons_root.rglob("*.polished.consensus.fasta")):
        sample = fasta.parent.name
        base = fasta.name
        suffix = ".polished.consensus.fasta"
        acc = ""
        if base.endswith(suffix):
            core = base[:-len(suffix)]
            prefix = f"{sample}."
            if core.startswith(prefix):
                acc = core[len(prefix):]
            else:
                parts = core.split(".", 1)
                acc = parts[1] if len(parts) == 2 else core

        total, called, n_bases, gap_bases, acgt = parse_fasta_counts(fasta)
        breadth = safe_pct(called, total)
        strict_breadth = safe_pct(acgt, total)
        rec = {
            "sample": sample,
            "accession": acc,
            "consensus_breadth_pct": breadth,
            "consensus_acgt_breadth_pct": strict_breadth,
        }
        rows.append(rec)
        by_pair[(sample, acc)] = rec

    qc_reads_by_sample = {}
if readcount_csv.exists():
    with readcount_csv.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            sample = (row.get("sample") or "").strip()
            trimmed = (row.get("trimmed_reads") or "").strip()
            if sample:
                qc_reads_by_sample[sample] = trimmed

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
    consensus_header = "Consensus Breadth (%)"

    out_fields = list(in_fields)

    # Remove any old technical consensus header so we can re-add consistently.
    out_fields = [f for f in out_fields if f != "consensus_breadth_pct"]

    # Enforce QC column before mapped reads.
    if qc_header in out_fields:
        out_fields.remove(qc_header)
    if "Mapped Reads (#)" in out_fields:
        mapped_idx = out_fields.index("Mapped Reads (#)")
        out_fields.insert(mapped_idx, qc_header)
    else:
        out_fields.append(qc_header)

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
            row[qc_header] = qc_reads_by_sample.get(sample, "")
            row[consensus_header] = "" if hit is None else hit["consensus_breadth_pct"]
    else:
        for row in data:
            row[qc_header] = ""
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

    inj_script = '''
<script>
(function(){
  const consensusMap = __CONS_MAP_JSON__;

  function getVal(sampleName, accession){
    if (!consensusMap || !consensusMap[sampleName]) return null;
    const v = consensusMap[sampleName][accession];
    if (v === undefined || v === null || v === '') return null;
    return Number(v);
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
    //    Wanted compact view:
    //      Accession | Organism | Viral species | Clade | Segment | Breadth | Consensus Breadth
    //    So hide: Mapped reads (#), Identity, Depth, NOGR, Z-score.
    const hideMainSorts = ['mapped_reads','identity','depth','nogr_regions','zscore'];
    hideMainSorts.forEach(s => {
      const th = mainRow.querySelector('th[data-sort="' + s + '"]');
      if (th) th.style.display = visible ? 'none' : '';
    });

    const hidePlaceholders = new Set(['Min Reads','Min ID','Min Depth','Min NOGR','Min Z-score']);
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
          tds[9].style.display = visible ? 'none' : '';
          tds[10].style.display = visible ? 'none' : '';
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

