process VIRASIGN_CLASSIFICATION {
    tag "Viral classification"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'library://jansendaan94_v2/metatropics/virasign:latest':
        'daanjansen94/virasign:latest' }"

    // Bind outdir/db_dir so results persist outside the container.
    def outAbs = file(params.outdir).toAbsolutePath().toString()
    def pipelineRoot = new File("${projectDir}").parentFile.absolutePath
    def db_dir = params.virasign_db_dir ?: "${pipelineRoot}/Databases"
    def dbAbs = file(db_dir).toAbsolutePath().toString()

    // Ensure any user-provided control/blind files are visible in-container.
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
    // One path per sample; used as a light barrier for METATROPICS_SUMMARY (avoid collect() on all publish/**).
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
    // Output isolation per database choice (avoid mixing when using -resume).
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

    # If -resume cached the barrier but DB files were deleted, retry prepare-db.
    if [ "${effectiveDbArg.toLowerCase()}" = "rvdb" ] && [ "${resolvedDbArg}" = "${effectiveDbArg}" ]; then
      virasign --prepare-db --db-dir "${db_dir}" -d 'RVDB' \
        ${params.virasign_rvdb_version != null ? "--rvdb-version ${params.virasign_rvdb_version}" : ""} \
        ${params.virasign_accessions?.toString()?.trim() ? "-a '${params.virasign_accessions}'" : ""} || true
    fi
    if [ "${effectiveDbArg.toLowerCase()}" = "refseq" ] && [ "${resolvedDbArg}" = "${effectiveDbArg}" ]; then
      virasign --prepare-db --db-dir "${db_dir}" -d 'RefSeq' \
        ${params.virasign_accessions?.toString()?.trim() ? "-a '${params.virasign_accessions}'" : ""} || true
    fi

    mkdir -p virasign_in

    # Stage input FASTQs into virasign_in/ (accept file or dir). Prefer canonical sample ID.
    stage_one () {
      local f="\$1"
      [ -e "\$f" ] || return 0

      local base="\${f##*/}"
      base="\${base%.fastq.gz}"
      base="\${base%.fq.gz}"
      base="\${base%_other}"

      local out="${meta.id}.fastq.gz"
      # If multiple files are provided, avoid name collisions.
      if [ -e "virasign_in/\$out" ]; then
        out="\${base}.fastq.gz"
      fi
      # Use absolute symlinks (workdir inputs may be symlinks themselves).
      target=\$(readlink -f "\$f" 2>/dev/null || realpath "\$f" 2>/dev/null || echo "\$f")
      ln -sfn "\$target" "virasign_in/\$out"
    }

    for f in ${virasign_input_fastqs}; do
      [ -e "\$f" ] || continue
      if [ -d "\$f" ]; then
        # Link all FASTQs in the directory.
        while IFS= read -r -d '' fq; do
          stage_one "\$fq"
        done < <(find "\$f" -maxdepth 1 -type f \\( -name '*.fastq' -o -name '*.fq' -o -name '*.fastq.gz' -o -name '*.fq.gz' \\) -print0)
      else
        stage_one "\$f"
      fi
    done

    n=\$(find virasign_in -maxdepth 1 -type l \\( -name '*.fastq' -o -name '*.fq' -o -name '*.fastq.gz' -o -name '*.fq.gz' \\) | wc -l)
    if [ "\$n" -eq 0 ]; then
      echo "ERROR: no FASTQs staged into virasign_in (input resolved to a directory with no FASTQ files, or empty optional upstream output)." >&2
      ls -la virasign_in >&2 || true
      exit 1
    fi
    ${cmd}

    # Create per-virus depth mask BEDs alongside Virasign BAMs.
    # These are later used to N-mask consensus sequences (depth < params.depth).
    # Virasign's container includes samtools.
    MIN_DEPTH="${params.depth ?: 25}"
    if command -v samtools >/dev/null 2>&1; then
      while IFS= read -r -d '' bam; do
        bai="\${bam}.bai"
        if [ ! -e "\$bai" ] && [ -e "\${bam%.bam}.bam.bai" ]; then
          bai="\${bam%.bam}.bam.bai"
        fi
        samtools index "\$bam" >/dev/null 2>&1 || true
        samtools depth -aa -d 0 "\$bam" \
          | awk -v min="\$MIN_DEPTH" 'BEGIN{OFS="\\t"} { if(\$3 < min) print \$1, \$2-1, \$2 }' \
          > "\${bam%.bam}.depth_lt_min.bed"
      done < <(find publish -type f -name '*.bam' -print0)
    else
      echo "WARNING: samtools not found in Virasign container; depth mask BEDs will not be generated." >&2
    fi

    # Don't propagate per-run virasign log into shared results.
    rm -f publish/.virasign.log || true

    # Copy publish/ into the shared results dir (parallel-safe with flock if available).
    mkdir -p "${shared}"
    if command -v flock >/dev/null 2>&1; then
      (
        flock -x 9
        for d in "\${PWD}/publish/"*; do
          [ -e "\$d" ] || continue
          bn="\$(basename \"\$d\")"
          echo "[virasign] Copying \${bn} -> ${shared}/" >&2
          cp -aL "\$d" "${shared}/"
          [ -e "${shared}/\${bn}" ] || { echo "[virasign] ERROR: destination missing after copy: ${shared}/\${bn}" >&2; exit 1; }
          # Do not keep internal mask artifacts in the published results tree.
          # They are only needed downstream during the current workflow execution.
          find "${shared}/\${bn}" -type f -name '*.depth_lt_min.bed' -delete 2>/dev/null || true
        done
      ) 9>"${shared}/.copy.lock"
    else
      for d in "\${PWD}/publish/"*; do
        [ -e "\$d" ] || continue
        bn="\$(basename \"\$d\")"
        echo "[virasign] Copying \${bn} -> ${shared}/" >&2
        cp -aL "\$d" "${shared}/"
        [ -e "${shared}/\${bn}" ] || { echo "[virasign] ERROR: destination missing after copy: ${shared}/\${bn}" >&2; exit 1; }
        find "${shared}/\${bn}" -type f -name '*.depth_lt_min.bed' -delete 2>/dev/null || true
      done
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        virasign: \$(virasign --version 2>&1 | head -n1 || echo 'unknown')
    END_VERSIONS
    """
}

