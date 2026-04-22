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
    path "results_summary_*.html", emit: html
    path "results_summary_*.csv",  emit: csv
    path "versions.yml",           emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def rawDbArg = params.virasign_database?.toString()?.trim()
    def effectiveDbArg = rawDbArg ?: 'RVDB'
    def virasignDbLabel = effectiveDbArg.replaceAll(/[^A-Za-z0-9._-]+/, '_')
    def outroot = file("${params.outdir}/Classification/virasign/${virasignDbLabel}").toAbsolutePath().toString()
    def zscore = params.virasign_zscore ?: 'true'
    def zscore_controls = params.virasign_zscore_controls?.toString()?.trim()
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
      ts=\$(date +%Y%m%d_%H%M%S)
      csv="results_summary_\${ts}.csv"
      html="results_summary_\${ts}.html"
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
    cp -f "${outroot}"/results_summary_*.html .
    cp -f "${outroot}"/results_summary_*.csv .
    # Keep the shared Classification/virasign tree focused on per-sample outputs only.
    rm -f "${outroot}"/results_summary_*.html "${outroot}"/results_summary_*.csv "${outroot}"/.virasign.log || true

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        virasign: \$(virasign --version 2>&1 | head -n1 || echo 'unknown')
    END_VERSIONS
    """
}

