process DORADO_DEMULTIPLEXING {
    label 'process_medium'

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

    # Get list of barcode files
    ls demultiplexed/*_barcode*.fastq > list.txt

    # Extract barcode names and move files
    while read file; do
        barcode=\$(echo \$file | grep -o 'barcode[0-9]*')
        mv "\$file" "\${barcode}.fastq"
    done < list.txt

    # Move unclassified file
    mv demultiplexed/*_unclassified.fastq unclassified.fastq

    # Get version
    VERSION=\$(dorado --version 2>&1 | tail -n 1)

    # Create versions file
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dorado: \$VERSION
    END_VERSIONS
    """
}