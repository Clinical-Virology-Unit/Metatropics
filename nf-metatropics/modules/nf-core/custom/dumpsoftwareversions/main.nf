process CUSTOM_DUMPSOFTWAREVERSIONS {
    label 'process_single'

    conda "conda-forge::python=3.12 conda-forge::pyyaml=6.0.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://python:3.12-slim-bookworm' :
        'python:3.12-slim-bookworm' }"
    beforeScript 'pip install --root-user-action=ignore --disable-pip-version-check -q pyyaml==6.0.2'

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
