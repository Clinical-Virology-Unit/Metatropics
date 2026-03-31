process METAMAPS_MAP {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::metamaps=0.1.98102e9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'library://daanjansen94/metatropics/metamaps:v0.1':
        'daanjansen94/metamaps:v0.1' }"

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("*_results"), emit: metaclass
    tuple val(meta), path("*_results.meta"), emit: otherclassmeta
    tuple val(meta), path("*_results.meta.unmappedReadsLengths"), emit: metalength
    tuple val(meta), path("*_results.parameters"), emit: metaparameters

    // TODO nf-core: List additional required output channels/values here
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    metamaps mapDirectly --all --pi 95 -w 100 $args -q $input -o ${prefix}_classification_results --maxmemory 12

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        metamaps_maps: \$(echo \$(metamaps --help) | grep MetaMaps | perl -p -e 's/MetaMaps v |Simul.+//g' )
    END_VERSIONS
    """
    }
