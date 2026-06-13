process VIRASIGN_CLASSIFICATION {
    tag "Viral classification"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'library://jansendaan94_v2/metatropics/virasign:latest':
        'daanjansen94/virasign:latest' }"

    // Bind db_dir and optional control files; virasign outputs reach outdir via publishDir (host-side copy).
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
        containerOptions "-v ${dbAbs}:${dbAbs}${extra}"
    }
    if (workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer') {
        def extra = extraBindDirs.collect { " --bind ${it}:${it}" }.join('')
        containerOptions "--bind ${dbAbs}:${dbAbs}${extra}"
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
    path "hits.tsv",                             emit: hits_tsv
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

    # Create per-virus mask BEDs alongside Virasign BAMs.
    # These are later used to N-mask consensus sequences.
    #
    # IMPORTANT: We enforce:
    # - minimum depth (`params.depth`, default 25)
    # - minimum base quality and mapping quality matching Clair3 uniform recount
    #   (`params.clair3_min_bq` / `params.clair3_min_mq`)
    #
    # This keeps site masking conceptually aligned with variant calling thresholds
    # (Clair3 uses --qual / --min_mq during calling; we at least mirror base-Q here).
    MIN_DEPTH="${params.depth ?: 25}"
    MIN_BQ="${params.clair3_min_bq ?: 15}"
    MIN_MQ="${params.clair3_min_mq ?: 15}"
    if command -v samtools >/dev/null 2>&1; then
      while IFS= read -r -d '' bam; do
        bai="\${bam}.bai"
        if [ ! -e "\$bai" ] && [ -e "\${bam%.bam}.bam.bai" ]; then
          bai="\${bam%.bam}.bam.bai"
        fi
        samtools index "\$bam" >/dev/null 2>&1 || true
        # Compute depth after filtering low-quality bases:
        # -Q MIN_BQ filters bases by Phred base quality
        # mpileup depth is column 4 (coverage) in samtools output
        ref="\${bam%.bam}.fasta"
        if [ ! -e "\$ref" ]; then
          # Fallback: keep old behaviour (depth-only) if reference isn't present for mpileup.
          samtools depth -aa -d 0 "\$bam" \
            | awk -v min="\$MIN_DEPTH" 'BEGIN{OFS="\\t"} { if(\$3 < min) print \$1, \$2-1, \$2 }' \
            > "\${bam%.bam}.bed"
        else
          samtools mpileup -aa -d 0 -Q "\$MIN_BQ" -q "\$MIN_MQ" -f "\$ref" "\$bam" 2>/dev/null \
            | awk -v min="\$MIN_DEPTH" 'BEGIN{OFS="\\t"} { if(\$4 < min) print \$1, \$2-1, \$2 }' \
            > "\${bam%.bam}.bed"
        fi
      done < <(find publish -type f -name '*.bam' -print0)
    else
      echo "WARNING: samtools not found in Virasign container; depth mask BEDs will not be generated." >&2
    fi

    # One TSV row per confident hit (same folder layout as Virasign: publish/<sample>/<acc>/).
    python3 <<PY
    import json, os, re, sys

    publish_dir = "publish"
    sample_id = "${meta.id}"
    sample_dir = os.path.join(publish_dir, sample_id)
    candidates = []
    if os.path.isdir(sample_dir):
        candidates.append(os.path.join(sample_dir, "%s_final_selected_references.json" % sample_id))
        for fn in os.listdir(sample_dir):
            if fn.endswith("_final_selected_references.json"):
                candidates.append(os.path.join(sample_dir, fn))
    candidates = [c for c in candidates if os.path.exists(c)]
    if not candidates:
        open("hits.tsv", "w").close()
        sys.exit(0)
    final_json = candidates[0]

    with open(final_json) as fh:
        hits = json.load(fh)

    def slug(s):
        s = (s or "").strip()
        if not s:
            return ""
        s = re.sub(r"[^A-Za-z0-9._-]+", "_", s)
        s = s.lstrip("_").rstrip("_")
        s = re.sub(r"_+", "_", s)
        return s

    out = open("hits.tsv", "w")
    for hit in hits:
        acc = str(hit.get("accession", "") or "").strip()
        if not acc:
            continue
        raw_sp = (hit.get("organism") or hit.get("viral_species") or "").strip()
        if (not raw_sp) and hit.get("description"):
            raw_sp = str(hit["description"]).strip()[:120]
        sp_slug = slug(raw_sp)
        virus_slug = "%s_%s" % (acc, sp_slug) if sp_slug else acc
        ref_fasta = os.path.join(os.path.dirname(final_json), acc, "%s.fasta" % acc)
        if not os.path.exists(ref_fasta):
            continue
        out.write("\\t".join([sample_id, acc, virus_slug, sp_slug, os.path.realpath(ref_fasta)]) + "\\n")
    out.close()
    PY

    # Don't propagate per-run virasign log into shared results.
    rm -f publish/.virasign.log || true

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        virasign: \$(virasign --version 2>&1 | head -n1 || echo 'unknown')
    END_VERSIONS
    """
}

