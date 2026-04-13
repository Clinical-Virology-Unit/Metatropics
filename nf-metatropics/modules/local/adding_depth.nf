process ADDING_DEPTH {
    tag "${meta.id}.${meta.virus}"

    input:
    tuple val(meta), path(depth), path(consensus), path(bamstats), path(report)

    output:
    tuple val(meta), path("*.sdepth.tsv"), emit : repdepth

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}.${meta.virus}"
    """
    adding_bcfdepth_V2.pl $depth $report ${prefix}.sdepth.tsv $consensus $bamstats
    """
}
