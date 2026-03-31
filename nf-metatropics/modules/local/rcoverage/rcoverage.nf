process RCOVERAGE {
    tag "rcoverage"
    label 'process_high' 

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'library://daanjansen94/metatropics/rcoverage:v4.0.2':
        'rocker/tidyverse:latest' }"

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
