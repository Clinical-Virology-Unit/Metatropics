process DORADO_DEMULTIPLEXING {
    tag "Demultiplexing"
    label 'process_gpu'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'library://daanjansen94/metatropics/dorado:0.9.0' :
        'nanoporetech/dorado:sha4644018526d3644d92a9e680ab8f2d1eeff2e272' }"

    // Set container options for Docker and Singularity
    containerOptions {
        if (workflow.containerEngine == 'singularity') {
            return '--nv --no-home'
        } else if (workflow.containerEngine == 'docker') {
            return '--gpus all --rm --init'
        } else {
            return null
        }
    }

    input:
    path reads

    output:
    path "*.fastq", emit: demultiplexed_fastq
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    dorado demux --kit-name ${params.kit_name} --emit-fastq --barcode-both-ends --output-dir demultiplexed $reads

    # Collect barcode outputs.
    ls demultiplexed/*_barcode*.fastq > list.txt

    # Rename by barcode.
    while read file; do
        barcode=\$(echo \$file | grep -o 'barcode[0-9]*')
        mv "\$file" "\${barcode}.fastq"
    done < list.txt

    # Rename unclassified.
    mv demultiplexed/*_unclassified.fastq unclassified.fastq

    # Version.
    VERSION=\$(dorado --version 2>&1 | tail -n 1)

    # Versions file.
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dorado: \$VERSION
    END_VERSIONS
    """
}
