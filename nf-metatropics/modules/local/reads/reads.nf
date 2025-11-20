process ReadCount {
    publishDir "${params.outdir}", mode: 'copy', overwrite: true, saveAs: { filename ->
        if (filename.endsWith('read_counts.csv')) return filename
        else if (filename.endsWith('read_distribution.pdf')) return filename
        else null
    }
    label 'process_medium'
    tag "ReadCount"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-tidyverse:1.2.1':
        'rocker/tidyverse:latest' }"

    if( workflow.containerEngine == 'docker' ) {
        def outPath = file(params.outdir).toAbsolutePath().toString()
        containerOptions "-v ${outPath}:${outPath} -u \$(id -u):\$(id -g)"
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

    # Copy raw reads from 'fix' folder when available
    if [ -d "${outdir}/fix" ]; then
        find ${outdir}/fix -name "*.fastq.gz" -type f -exec cp {} read_count/ \\;
    else
        echo "Directory ${outdir}/fix does not exist, skipping raw read copy"
    fi

    # Copy trimmed reads from 'fastp' folder when available
    if [ -d "${outdir}/fastp" ]; then
        find ${outdir}/fastp -name "*.fastq.gz" -type f -exec cp {} read_count/ \\;
    else
        echo "Directory ${outdir}/fastp does not exist, skipping trimmed read copy"
    fi

    # Copy human-depleted reads from 'nohuman' folder when available
    if [[ "\$HOST_STATUS" == "human_only" || "\$HOST_STATUS" == "both" ]]; then
        if [ -d "${outdir}/nohuman" ]; then
            find ${outdir}/nohuman -name '*other.fastq.gz' -type f -exec cp {} read_count/nohuman/ \\;
        else
            echo "Directory ${outdir}/nohuman does not exist, skipping human-depleted copy"
        fi
    else
        echo "Human host depletion not enabled; skipping human-depleted copy"
    fi

    # Copy host-depleted reads from 'nohost' folder if it exists
    if [[ "\$HOST_STATUS" == "other_only" || "\$HOST_STATUS" == "both" ]]; then
        if [ -d "${outdir}/nohost" ]; then
            find ${outdir}/nohost -name '*other.fastq.gz' -type f -exec cp {} read_count/nohost/ \\;
        else
            echo "Directory ${outdir}/nohost does not exist, skipping host-depleted copy"
        fi
    else
        echo "Additional host depletion not enabled; skipping host-depleted copy"
    fi

    # Process viral reads from 'metamaps' folder
    echo "Sample,ViralReads" > read_count/viral_read_counts.csv
    if [ -d "${outdir}/metamaps" ]; then
        found_files=false
        for file in ${outdir}/metamaps/*_classification_results.meta; do
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
        echo "Directory ${outdir}/metamaps does not exist, skipping viral read count extraction"
    fi

    ReadCount.R ${params.outdir}/read_count/ ${host_genome_status}
    """
}
