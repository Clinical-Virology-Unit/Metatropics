process METAMAPS_CLASSIFY {
    tag "$meta.id"
    label 'process_high'
    conda "bioconda::metamaps=0.1.98102e9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'library://daanjansen94/metatropics/metamaps:v0.1':
        'daanjansen94/metamaps:v0.1' }"

    input:
    tuple val(meta), path(input), path(metamap), path(unmapped), path(parametersmeta)

    output:
    tuple val(meta), path("*_results.unique_virus"), emit: classem
    tuple val(meta), path("*_results.EM"), emit: classem_original
    tuple val(meta), path("*.EM.reads2Taxon.krona"), emit: classkrona
    tuple val(meta), path("*.EM.lengthAndIdentitiesPerMappingUnit"), emit: classlength
    tuple val(meta), path("*.EM.WIMP"), emit: classWIMP
    tuple val(meta), path("*.EM.contigCoverage"), emit: classcov
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Run MetaMaps classify
    # Note: metamaps may crash with exit 134 (assertion failure) after completing EM algorithm
    # The .EM file is created before the crash, so we can use it even if the command fails
    set +e
    metamaps classify --mappings ${prefix}_classification_results $args
    classify_exit=\$?
    set -e
    
    # Check if the essential output file was created (even if command failed)
    if [ ! -f "${prefix}_classification_results.EM" ]; then
        echo "ERROR: ${prefix}_classification_results.EM was not created. Metamaps classify failed before producing output."
        exit 1
    fi
    
    # If metamaps crashed with exit 134 but .EM file exists, log a warning but continue
    # This is a known bug where metamaps crashes after completing the EM algorithm
    if [ \$classify_exit -eq 134 ]; then
        echo "WARNING: Metamaps classify crashed with assertion failure (exit 134), but .EM file exists."
        echo "WARNING: Using existing .EM file to continue processing. This is a known metamaps bug."
        echo "WARNING: The classification results should still be valid."
        # Exit with 0 to indicate success since we have the output we need
        # (we'll re-enable strict error checking after this)
    elif [ \$classify_exit -ne 0 ]; then
        echo "ERROR: Metamaps classify failed with exit code \$classify_exit"
        exit \$classify_exit
    fi

    # Filter for unambiguous mappings
    awk '
    {
        read_id = \$1
        alignment_score = \$10
        mapping_quality = \$NF
        combined_score = alignment_score * mapping_quality
        if (!(read_id in best_score) || combined_score > best_score[read_id]) {
            best_score[read_id] = combined_score
            best_line[read_id] = \$0
        }
    }
    END {
        for (read_id in best_line) {
            print best_line[read_id]
        }
    }' ${prefix}_classification_results.EM > ${prefix}_classification_results.unique_virus

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        metamaps_classify: \$(echo \$(metamaps --help) | grep MetaMaps | perl -p -e 's/MetaMaps v |Simul.+//g' )
    END_VERSIONS
    """
}
