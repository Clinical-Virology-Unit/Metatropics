process ReadCount {
    label 'process_medium'
    tag "Read distribution summary"
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

    # Copy raw reads from Reads/raw when available
    if [ -d "${outdir}/Reads/raw" ]; then
        find ${outdir}/Reads/raw -name "*.fastq.gz" -type f -exec cp {} read_count/ \\;
    else
        echo "Directory ${outdir}/Reads/raw does not exist, skipping raw read copy"
    fi

    # Copy trimmed reads from Reads/fastplong when available
    if [ -d "${outdir}/Reads/fastplong" ]; then
        find ${outdir}/Reads/fastplong -name "*.fastq.gz" -type f -exec cp {} read_count/ \\;
    else
        echo "Directory ${outdir}/Reads/fastplong does not exist, skipping trimmed read copy"
    fi

    # Copy human-depleted reads (new or legacy naming).
    if [[ "\$HOST_STATUS" == "human_only" || "\$HOST_STATUS" == "both" ]]; then
        if [ -d "${outdir}/Reads/nohuman" ]; then
            find ${outdir}/Reads/nohuman \\( -name '*_human_depleted.fastq.gz' -o -name '*_other.fastq.gz' \\) -type f -exec cp {} read_count/nohuman/ \\;
        else
            echo "Directory ${outdir}/Reads/nohuman does not exist, skipping human-depleted copy"
        fi
    else
        echo "Human host depletion not enabled; skipping human-depleted copy"
    fi

    # Copy host-depleted reads (new or legacy naming).
    if [[ "\$HOST_STATUS" == "other_only" || "\$HOST_STATUS" == "both" ]]; then
        if [ -d "${outdir}/Reads/nohost" ]; then
            find ${outdir}/Reads/nohost \\( -name '*_host_depleted.fastq.gz' -o -name '*_other.fastq.gz' \\) -type f -exec cp {} read_count/nohost/ \\;
        else
            echo "Directory ${outdir}/Reads/nohost does not exist, skipping host-depleted copy"
        fi
    else
        echo "Additional host depletion not enabled; skipping host-depleted copy"
    fi

    # Generate read count outputs.
    # Install to a writable per-task prefix so this works in read-only containers (Apptainer/Singularity).
    export PYTHONUSERBASE="\$PWD/.pyuserbase"
    mkdir -p "\$PYTHONUSERBASE"
    python -m pip install --no-cache-dir --user pandas plotly matplotlib >/dev/null
    python ${projectDir}/bin/readcount.py \
      --outdir "${outdir}" \
      --host-status "${host_genome_status}" \
      --workdir "."

    # Clean up staged FASTQs (counting only).
    find read_count -type f -name '*.fastq.gz' -delete || true
    """
}
