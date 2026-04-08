process FIX_NAMES {
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'library://daanjansen94/metatropics/bbmap:38.86':
        'daanjansen94/bbmap:38.86' }"

    tag{sample}

    input:
    tuple val(meta), val(sample), path(reads)

    output:
    tuple val(sample), path("*.fastq.gz"), emit : fqreads

    script:
    """
    cat $reads > ${sample}_fixed.fastq
    reformat.sh in=${sample}_fixed.fastq out=${sample}_fixed.fastq.gz qin=33 ignorebadquality overwrite=t
    """
}
