process NANOPLOT {
    tag "Read QC (NanoPlot)"
    label 'process_single'

    conda "bioconda::nanoplot=1.41.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'library://daanjansen94/metatropics/nanoplot:v1.41.0' :
        'daanjansen94/nanoplot:v1.41.0' }"

    input:
    tuple val(meta), path(ontfile)

    output:
    tuple val(meta), path("*.html")                , emit: html
    tuple val(meta), path("*.png") , optional: true, emit: png
    tuple val(meta), path("*.txt")                 , emit: txt
    tuple val(meta), path("*.log")                 , emit: log
    tuple val(meta), path("*.total_reads")         , emit: totalreads
    path  "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def ont = ontfile.toString()
    def input_file = (ont.endsWith('.fastq.gz') || ont.endsWith('.fq.gz')) ? "--fastq ${ontfile}" :
        (ont.endsWith('.fastq') || ont.endsWith('.fq')) ? "--fastq ${ontfile}" :
        (ont.endsWith('.txt')) ? "--summary ${ontfile}" : ''
    """
    if gzip -t $ontfile 2>/dev/null; then
        gzip -cd -- $ontfile | wc -l | awk '{x=\$1/4; print x}' > ${meta.id}_classification_results.total_reads
    else
        wc -l < $ontfile | awk '{x=\$1/4; print x}' > ${meta.id}_classification_results.total_reads
    fi
    NanoPlot \\
        $args \\
        -t $task.cpus \\
        -p $meta.id \\
        $input_file
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanoplot: \$(echo \$(NanoPlot --version 2>&1) | sed 's/^.*NanoPlot //; s/ .*\$//')
    END_VERSIONS
    """
}
