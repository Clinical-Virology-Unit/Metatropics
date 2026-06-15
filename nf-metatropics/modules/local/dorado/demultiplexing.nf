process DORADO_DEMULTIPLEXING {
    tag "Demultiplexing"
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

    # Collect barcode FASTQs (Dorado 0.x flat layout and 2.x nested fastq_pass layout).
    find demultiplexed -type f \\( -path '*/fastq_pass/barcode*/*.fastq' -o -name '*_barcode*.fastq' \\) ! -path '*/unclassified/*' | sort -u > list.txt

    if [[ ! -s list.txt ]]; then
        echo "ERROR: No demultiplexed barcode FASTQ files found under demultiplexed/" >&2
        find demultiplexed -type f | head -20 >&2 || true
        exit 1
    fi

    # Rename by barcode.
    while IFS= read -r file; do
        barcode=\$(echo "\$file" | grep -oE 'barcode[0-9]+' | head -1)
        if [[ -n "\$barcode" ]]; then
            mv "\$file" "\${barcode}.fastq"
        fi
    done < list.txt

    # Rename unclassified (if present).
    unclassified=\$(find demultiplexed -type f \\( -path '*/fastq_pass/unclassified/*.fastq' -o -name '*_unclassified.fastq' \\) | head -1 || true)
    if [[ -n "\$unclassified" ]]; then
        mv "\$unclassified" unclassified.fastq
    fi

    # Version.
    VERSION=\$(dorado --version 2>&1 | tail -n 1)

    # Versions file.
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dorado: \$VERSION
    END_VERSIONS
    """
}
