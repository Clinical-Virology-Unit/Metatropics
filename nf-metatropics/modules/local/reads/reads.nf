process ReadCount {
    label 'process_medium'
    tag "ReadCount"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'library://jansendaan94_v2/metatropics/rocker_tidyverse:latest':
        'daanjansen94/rocker-tidyverse:latest' }"

    def outPath = file(params.outdir).toAbsolutePath().toString()
    if( workflow.containerEngine == 'docker' ) {
        containerOptions "-v ${outPath}:${outPath} -u \$(id -u):\$(id -g)"
    }
    if( workflow.containerEngine == 'singularity' ) {
        containerOptions "--bind ${outPath}:${outPath}"
    }

    input:
    val outdir
    path medaka_files
    val host_genome_status

    output:
    path "read_count/*.fastq.gz", emit: read_count_fastq_root, optional: true
    path "read_count/**/*.fastq.gz", emit: read_count_fastq_nested, optional: true
    path "read_count/read_counts.csv", emit: read_counts_csv
    path "read_count/read_distribution.pdf", emit: read_distribution_pdf

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

    # Copy trimmed reads from Reads/fastp when available
    if [ -d "${outdir}/Reads/fastp" ]; then
        find ${outdir}/Reads/fastp -name "*.fastq.gz" -type f -exec cp {} read_count/ \\;
    else
        echo "Directory ${outdir}/Reads/fastp does not exist, skipping trimmed read copy"
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

    # Process viral reads from Classification/metamaps
    echo "Sample,ViralReads" > read_count/viral_read_counts.csv
    if [ -d "${outdir}/Classification/metamaps" ]; then
        found_files=false
        for file in ${outdir}/Classification/metamaps/*_classification_results.meta; do
            if [ -f "\$file" ]; then
                found_files=true
                sample_name=\$(basename "\$file" _classification_results.meta)
                viral_reads=\$(grep "ReadsMapped" "\$file" | awk '{print \$2}')
                echo "\${sample_name},\${viral_reads}" >> read_count/viral_read_counts.csv
            fi
        done
        if [ "\$found_files" = false ]; then
            echo "No matching files found in metamaps folder"
        fi
    else
        echo "Directory ${outdir}/Classification/metamaps does not exist, skipping viral read count extraction"
    fi

    ReadCount.R ${params.outdir}/Summary/read_count/ ${host_genome_status}
    """
}
