process R_METAPLOT {
    tag "$meta.id"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'library://daanjansen94/metatropics/nf:r_plots_v4.2.2':
        'daanjansen94/nf_r_plots:v4.2.2' }"

    input:
    tuple val(meta), path(classification_results), path(length_and_identities), path(contig_coverage), path(total_reads), path(cleanup_done)

    output:
    tuple val(meta), path("*.pdf"), emit: plotpdf
    tuple val(meta), path("*.tsv"), emit: reporttsv
    tuple val(meta), path("*.txt"), emit: denovo
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Bundled script (overrides older copy that may exist in the container)
    Rscript ${projectDir}/bin/plotMappingSummary.R ${prefix}_classification_results $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R used to plots: \$(echo \$(R --version) | perl -p -e 's/R version //g' | perl -p -e 's/ .+//g' )
    END_VERSIONS
    """
}
