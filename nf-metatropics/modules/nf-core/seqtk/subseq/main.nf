process SEQTK_SUBSEQ {
    tag "$sequences"
    label 'process_single'

    conda "bioconda::seqtk=1.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'library://daanjansen94/metatropics/seqtk:v1.3' :
        'daanjansen94/seqtk:v1.3' }"

    input:
    tuple val(meta), path (sequences), path (filter_list)
    //path filter_list

    output:
    tuple val(meta), path("*.gz") , emit: sequences
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: ''
    def ext = "fa"
    if ("$sequences" ==~ /.+\.fq|.+\.fq.gz|.+\.fastq|.+\.fastq.gz/) {
        ext = "fq"
    }
    """
    seqtk \\
        subseq \\
        $args \\
        $sequences \\
        $filter_list | \\
        gzip --no-name > ${sequences}${prefix}.${ext}.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqtk: \$(echo \$(seqtk 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """
}
