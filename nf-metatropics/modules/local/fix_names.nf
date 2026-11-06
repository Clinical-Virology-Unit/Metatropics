process FIX_NAMES {
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://daanjansen94/bbmap:38.86':
        'daanjansen94/bbmap:38.86' }"

    tag "Evaluate/fix format of raw reads"

    input:
    tuple val(meta), val(sample), path(reads)

    output:
    tuple val(sample), path("${sample}_fixed.fastq.gz"), emit : fqreads

    script:
    // barcode column may be one FASTQ file, or a folder of FASTQ shards
    """
    set -euo pipefail
    if [ -d "${reads}" ]; then
        shopt -s nullglob
        files=( "${reads}"/*.fastq.gz "${reads}"/*.fq.gz "${reads}"/*.fastq "${reads}"/*.fq )
        if [ \${#files[@]} -eq 0 ]; then
            echo "ERROR: no FASTQ files under ${reads}" >&2
            exit 1
        fi
        if [ \${#files[@]} -eq 1 ]; then
            in="\${files[0]}"
        else
            cat "\${files[@]}" > merged_input.fastq.gz
            in=merged_input.fastq.gz
        fi
    else
        in="${reads}"
    fi
    reformat.sh in="\$in" out=${sample}_fixed.fastq.gz qin=33 ignorebadquality overwrite=t
    rm -f merged_input.fastq.gz
    """
}
