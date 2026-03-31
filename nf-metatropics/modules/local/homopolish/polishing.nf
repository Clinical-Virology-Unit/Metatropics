process HOMOPOLISH_POLISHING {
    tag "${meta.id}.${meta.virus}"
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
    // TODO nf-core: Where possible, a command MUST be provided to obtain the version number of the software e.g. 1.10
    //               If the software is unable to output a version number on the command-line then it can be manually specified
    //               e.g. https://github.com/nf-core/modules/blob/master/modules/nf-core/homer/annotatepeaks/main.nf
    //               Each software used MUST provide the software name and version number in the YAML version file (versions.yml)
    // TODO nf-core: It MUST be possible to pass additional parameters to the tool as a command-line string via the "task.ext.args" directive
    // TODO nf-core: If the tool supports multi-threading then you MUST provide the appropriate parameter
    //               using the Nextflow "task" variable e.g. "--threads $task.cpus"
    // TODO nf-core: Please replace the example samtools command below with your module's command
    // TODO nf-core: Please indent the command appropriately (4 spaces!!) to help with readability ;)
    """
    homopolish polish -a $consensus -l $reffasta $args -o $prefix
    mv $prefix/* ${prefix}.polish.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        homopolish: \$(echo \$(homopolish --version) | perl -p -e 's/Homo.+: //g'  )
    END_VERSIONS
    """
}
