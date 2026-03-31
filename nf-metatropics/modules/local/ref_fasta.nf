process REF_FASTA {
    tag "$meta.id"

    conda "bioconda::metamaps=0.1.98102e9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'library://daanjansen94/metatropics/samtools:minimap2_v1.16_v2.28':
        'daanjansen94/minimap:v2.28' }"

    input:
    tuple val(meta), path(report), path(emreads), path(rawfastq)

    output:
    tuple val(meta), path("*.fasta"), emit : seqref
    tuple val(meta), path("*.reads"), emit : headereads
    tuple val(meta), path("*.fastq"), emit : allreads
    //tuple val(meta), stdout, emit : virusout

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    produce_fasta.pl $report $emreads $rawfastq $args
    """
}
