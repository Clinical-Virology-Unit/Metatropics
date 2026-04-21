process VIRASIGN_CLASSIFICATION {
    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'library://jansendaan94_v2/metatropics/virasign:latest':
        'daanjansen94/virasign:latest' }"

    // Bind host outdir so the process can merge results into a stable shared path.
    // Passing the shared path as a `path` input has intermittently serialized to `false` on some setups.
    def outAbs = file(params.outdir).toAbsolutePath().toString()
    def pipelineRoot = new File("${projectDir}").parentFile.absolutePath
    def db_dir = params.virasign_db_dir ?: "${pipelineRoot}/Databases"
    def dbAbs = file(db_dir).toAbsolutePath().toString()

    // If the user provides absolute paths for z-score controls / blind list, those files must
    // be visible inside the container as well. Bind their parent directories.
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
        containerOptions "-v ${outAbs}:${outAbs} -v ${dbAbs}:${dbAbs}${extra}"
    }
    if (workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer') {
        def extra = extraBindDirs.collect { " --bind ${it}:${it}" }.join('')
        containerOptions "--bind ${outAbs}:${outAbs} --bind ${dbAbs}:${dbAbs}${extra}"
    }

    input:
    path db_ready
    tuple val(meta), path(virasign_input_fastqs)

    output:
    path "publish/**",                          emit: results
    // One path per sample; used as a light barrier for VIRASIGN_SUMMARY (avoid collect() on all publish/**).
     path "publish/**/*_final_selected_references.json", emit: final_json, optional: true
    // If no confident hits, Virasign still writes the unfiltered JSON; treat that as a valid completion signal too.
    path "publish/**/*_unfiltered_all_references.json", emit: unfiltered_json, optional: true
    tuple val(meta), path("publish"),           emit: outdir
    path "versions.yml",                       emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Ensure host DB dir exists before the container writes into it.
    file(db_dir).mkdirs()
    def threads = params.virasign_threads ?: task.cpus
    def opt = []
    def add = { c, f -> if (c) opt << f }

    def rawDbArg = params.virasign_database?.toString()?.trim()
    def effectiveDbArg = rawDbArg ?: 'RVDB'
    // Database label for output isolation. This prevents mixing results across runs when
    // users change `virasign_database` but run with `-resume` (cached processes can
    // otherwise leave older files under a shared Classification/virasign root).
    def virasignDbLabel = effectiveDbArg.replaceAll(/[^A-Za-z0-9._-]+/, '_')
    def resolvedDbArg = effectiveDbArg
    if (effectiveDbArg) {
        def lower = effectiveDbArg.toLowerCase()
        if (lower == 'refseq') {
            def refseqFasta = file("${db_dir}/RefSeq/viral_refseq_complete.fna")
            if (refseqFasta.exists()) {
                resolvedDbArg = refseqFasta.toAbsolutePath().toString()
            }
        } else if (lower == 'rvdb') {
            def rvdbDir = file("${db_dir}/RVDB")
            if (rvdbDir.exists()) {
                def candidates = rvdbDir.listFiles()?.findAll { it.name ==~ /^RVDB.*_complete\.fasta$/ }
                if (candidates) {
                    resolvedDbArg = candidates.sort { it.name }.last().toAbsolutePath().toString()
                }
            }
        }
    }

    add(resolvedDbArg, "-d '${resolvedDbArg}'")
    add(params.virasign_rvdb_version != null, "--rvdb-version ${params.virasign_rvdb_version}")
    add(params.virasign_accessions?.toString()?.trim(), "-a '${params.virasign_accessions}'")
    add(params.virasign_ultrasensitive == true, '-u')
    add(params.virasign_min_identity != null, "--min_identity ${params.virasign_min_identity}")
    add(params.virasign_min_mapped_reads != null, "--min_mapped_reads ${params.virasign_min_mapped_reads}")
    add(params.virasign_coverage_depth != null, "--coverage_depth ${params.virasign_coverage_depth}")
    add(params.virasign_coverage_breadth != null, "--coverage_breadth ${params.virasign_coverage_breadth}")
    add(params.virasign_min_nogr != null, "--min-nogr ${params.virasign_min_nogr}")
    add(true, '--no-html')
    add(params.virasign_no_gzip_fastq == true, '--no-gzip-fastq')
    add(!!params.virasign_zscore, "--zscore ${params.virasign_zscore}")
    add(params.virasign_zscore_controls?.toString()?.trim(), "--zscore-controls '${params.virasign_zscore_controls}'")
    add(params.virasign_blind?.toString()?.trim(), "-b '${params.virasign_blind}'")

    def tail = task.ext.args?.toString()?.trim()
    def cmd = (['virasign', '-i', 'virasign_in', '-o', 'publish', "--db-dir", "${db_dir}", '-t', "${threads}"] + opt + (tail ? [tail] : [])).join(' ')
    def shared = file("${params.outdir}/Classification/virasign/${virasignDbLabel}").toAbsolutePath().toString()
    """
    # Barrier input from DB prep (ensures DB is prepared once).
    test -f "${db_ready}"

    mkdir -p virasign_in
    for f in *.fastq.gz *.fq.gz; do
      [ -e "\$f" ] || continue
      [ -d "\$f" ] && continue
      ln -sfn "\$PWD/\$f" "virasign_in/\$f"
    done
    n=\$(find virasign_in -mindepth 1 -maxdepth 1 | wc -l)
    if [ "\$n" -eq 0 ]; then
      echo "ERROR: no FASTQs staged into virasign_in (expected symlinks *.fastq.gz or *.fq.gz in work dir)." >&2
      ls -la >&2
      exit 1
    fi
    ${cmd}

    # Don't propagate the per-run virasign log into the shared results tree.
    rm -f publish/.virasign.log || true

    # Copy this sample's publish/ tree into the shared results dir.
    # Use flock for parallel-safety, but lock the shared directory itself (no extra .lock file created).
    mkdir -p "${shared}"
    if command -v flock >/dev/null 2>&1; then
      flock -x "${shared}" bash -lc "
        for d in \"\\\$PWD\"/publish/*; do
          [ -e \"\\\$d\" ] || continue
          cp -aL \"\\\$d\" \"${shared}/\"
        done
      "
    else
      for d in "\${PWD}/publish/"*; do
        [ -e "\$d" ] || continue
        cp -aL "\$d" "${shared}/"
      done
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        virasign: \$(virasign --version 2>&1 | head -n1 || echo 'unknown')
    END_VERSIONS
    """
}

