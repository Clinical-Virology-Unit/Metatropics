process CUSTOM_DUMPSOFTWAREVERSIONS {
    tag "Collect software versions"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'library://daanjansen94/metatropics/multiqc:v1.12' :
        'daanjansen94/multiqc:v1.12' }"

    input:
    path versions
    
    output:
    path "software_versions.yml"    , emit: yml
    path "software_versions_mqc.yml", emit: mqc_yml
    path "versions.yml"             , emit: versions
 
    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    template 'dumpsoftwareversions.py'
}
