process CUSTOM_DUMPSOFTWAREVERSIONS {
    label 'process_single'

    // Requires `pyyaml` which does not have a dedicated container but is in the MultiQC container
    conda "bioconda::multiqc=1.13"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'library://daanjansen94/metatropics/multiqc:v1.12' :
        'daanjansen94/multiqc:v1.12' }"

    input:
    path versions
    path rcoverage_done
    
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
