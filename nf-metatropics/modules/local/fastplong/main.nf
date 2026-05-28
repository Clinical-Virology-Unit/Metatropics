process FASTPLONG {
    tag "Read quality control"
    label 'process_medium'

    conda "bioconda::fastplong=0.4.1"
    // Singularity: Sylabs image lives under jansendaan94_v2 (same account as other metatropics SIFs you publish). Docker Hub: daanjansen94/fastplong
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'library://jansendaan94_v2/metatropics/fastplong:0.4.1' :
        'daanjansen94/fastplong:0.4.1' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*.fastp.fastq.gz'), emit: reads
    tuple val(meta), path('*.json'),        emit: json
    tuple val(meta), path('*.html'),        emit: html
    tuple val(meta), path('*.log'),         emit: log
    path "versions.yml",                    emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    if (!meta.single_end) {
        error 'FASTPLONG supports single-end (typical ONT) reads only. Paired-end samples are not supported by this module.'
    }
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reads_name = new File(reads.toString()).name
    def reads_stem = reads_name.endsWith('.gz') ? "${prefix}.fastq.gz" : "${prefix}.fastq"
    """
    set -euo pipefail
    [ ! -f ${reads_stem} ] && ln -sf $reads ${reads_stem}

    # Long-read QC (replaces fastp): same trim/QC intent; -A disables adapter trimming.
    fastplong \\
        -i ${reads_stem} \\
        -o ${prefix}.fastp.fastq.gz \\
        -z 4 \\
        -j ${prefix}.fastplong.json \\
        -h ${prefix}.fastplong.html \\
        -w $task.cpus \\
        -A \\
        $args \\
        2> ${prefix}.fastplong.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastplong: \$(fastplong --version 2>&1 | sed -e 's/fastplong //g')
    END_VERSIONS
    """

}
