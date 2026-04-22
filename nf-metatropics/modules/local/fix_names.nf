process FIX_NAMES {
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'library://daanjansen94/metatropics/bbmap:38.86':
        'daanjansen94/bbmap:38.86' }"

    tag "Evaluate/fix format of raw reads"

    input:
    tuple val(meta), val(sample), path(reads)

    output:
    tuple val(sample), path("*.fastq.gz"), emit : fqreads

    script:
    // reformat.sh reads gzip or plain FASTQ from path (avoid cat, which breaks .fastq.gz)
    """
    reformat.sh in=$reads out=${sample}_fixed.fastq.gz qin=33 ignorebadquality overwrite=t
    """
}
