process ReadCount {
    label 'process_medium'
    tag "ReadCount"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://python:3.11-bookworm' :
        'python:3.11-bookworm' }"

    def outPath = file(params.outdir).toAbsolutePath().toString()
    if( workflow.containerEngine == 'docker' ) {
        containerOptions "-v ${outPath}:${outPath} -u 0:0"
    }
    if( workflow.containerEngine == 'singularity' ) {
        containerOptions "--bind ${outPath}:${outPath}"
    }

    input:
    tuple val(outdir), val(_readcountBarrier), val(host_genome_status)

    output:
    path "read_count/read_counts.csv", emit: read_counts_csv
    path "read_count/read_distribution.pdf", emit: read_distribution_pdf
    path "read_count/read_distribution.html", emit: read_distribution_html

    script:
    """
    mkdir -p read_count read_count/nohuman read_count/nohost
    HOST_STATUS="${host_genome_status}"

    # Copy raw reads from Reads/fix when available
    if [ -d "${outdir}/Reads/fix" ]; then
        find ${outdir}/Reads/fix -name "*.fastq.gz" -type f -exec cp {} read_count/ \\;
    else
        echo "Directory ${outdir}/Reads/fix does not exist, skipping raw read copy"
    fi

    # Copy trimmed reads from Reads/fastplong when available
    if [ -d "${outdir}/Reads/fastplong" ]; then
        find ${outdir}/Reads/fastplong -name "*.fastq.gz" -type f -exec cp {} read_count/ \\;
    else
        echo "Directory ${outdir}/Reads/fastplong does not exist, skipping trimmed read copy"
    fi

    # Copy human-depleted reads from Reads/nohuman when available
    if [[ "\$HOST_STATUS" == "human_only" || "\$HOST_STATUS" == "both" ]]; then
        if [ -d "${outdir}/Reads/nohuman" ]; then
            find ${outdir}/Reads/nohuman -name '*other.fastq.gz' -type f -exec cp {} read_count/nohuman/ \\;
        else
            echo "Directory ${outdir}/Reads/nohuman does not exist, skipping human-depleted copy"
        fi
    else
        echo "Human host depletion not enabled; skipping human-depleted copy"
    fi

    # Copy host-depleted reads from Reads/nohost if it exists
    if [[ "\$HOST_STATUS" == "other_only" || "\$HOST_STATUS" == "both" ]]; then
        if [ -d "${outdir}/Reads/nohost" ]; then
            find ${outdir}/Reads/nohost -name '*other.fastq.gz' -type f -exec cp {} read_count/nohost/ \\;
        else
            echo "Directory ${outdir}/Reads/nohost does not exist, skipping host-depleted copy"
        fi
    else
        echo "Additional host depletion not enabled; skipping host-depleted copy"
    fi

    # Generate read_counts.csv + interactive HTML + PDF using the bundled dashboard script.
    python -m pip install --no-cache-dir pandas plotly matplotlib >/dev/null
    python ${projectDir}/bin/readcount.py \
      --outdir "${outdir}" \
      --host-status "${host_genome_status}" \
      --workdir "."

    # Staged FASTQs were only for counting; Summary publishes CSV/HTML/PDF only (see modules.config).
    find read_count -type f -name '*.fastq.gz' -delete || true
    """
}
