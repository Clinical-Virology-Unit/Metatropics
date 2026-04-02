process RCOVERAGE {
    tag "rcoverage"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'library://daanjansen94/metatropics/rcoverage:v4.0.2':
        'daanjansen94/rcoverage:v4.0.2' }"

    when:
    params.rcoverage_figure
    
    input:
    path coveragefiles

    output:
    path "coverage_distribution_group_*.pdf"

    script:
    """
    Coverage.R ${coveragefiles}
    """
}
