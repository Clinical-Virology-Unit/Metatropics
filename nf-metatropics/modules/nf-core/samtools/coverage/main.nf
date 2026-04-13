process SAMTOOLS_COVERAGE {
    tag "${meta.id}.${meta.virus}"
    label 'process_single'

    conda "bioconda::samtools=1.17"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'library://daanjansen94/metatropics/samtools:v1.17' :
        'daanjansen94/samtools:v1.17' }"

    input:
    tuple val(meta), path(input), path(input_index)

    output:
    tuple val(meta), path("*.txt"), emit: coverage
    tuple val(meta), path("*.bamstats"), emit: bamstats
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}.${meta.virus}"
    """
    samtools \\
        coverage \\
        $args \\
        -o ${prefix}.txt \\
        $input

    samtools view -F 2308 $input | awk '{
        nm = -1
        for(i=12; i<=NF; i++) {
            if(\$i ~ /^NM:i:/) {
                nm = substr(\$i, 6) + 0
                break
            }
        }
        if(nm < 0) next
        if(\$10 == "*") next
        qlen = length(\$10)
        if(qlen > 0) {
            identity = 1 - (nm / qlen)
            total_identity += identity
            total_length += qlen
            count++
        }
    }
    END {
        if(count > 0) {
            printf "%.6f\\t%.1f\\n", total_identity/count, total_length/count
        } else {
            printf "0\\t0\\n"
        }
    }' > ${prefix}.bamstats

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' )
    END_VERSIONS
    """
}
