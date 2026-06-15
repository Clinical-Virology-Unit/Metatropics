process DORADO_ONT {
    tag "Basecalling"
    label 'process_gpu'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://nanoporetech/dorado:latest' :
        'nanoporetech/dorado:latest' }"

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
    path input_dir

    output:
    path "calls.bam", emit: basecalling_ch
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    # Run dorado basecaller 
    dorado basecaller ${params.model} $input_dir \
        --no-trim \
        --batchsize 1024 \
        --chunksize 10000 \
        --overlap 500 \
        -x "auto" > calls.bam

    # Get version
    VERSION=\$(dorado --version 2>&1 | tail -n 1)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dorado: \$VERSION
    END_VERSIONS
    """
}

