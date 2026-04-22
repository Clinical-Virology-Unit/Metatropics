process HOMOPOLISH_POLISHING {
    tag "Consensus polishing"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'library://daanjansen94/metatropics/homopolish:v0.4.1':
        'daanjansen94/homopolish:v0.4.1' }"

    input:
    tuple val(meta), path(consensus), path(reffasta)

    output:
    tuple val(meta), path("*.polish.fasta"), emit: polishconsensus
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}.${meta.virus}"
    """
    homopolish polish -a $consensus -l $reffasta $args -o $prefix
    mv $prefix/* ${prefix}.polish.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        homopolish: \$(echo \$(homopolish --version) | perl -p -e 's/Homo.+: //g'  )
    END_VERSIONS
    """
}
