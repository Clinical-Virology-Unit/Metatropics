process SAMTOOLS_hFASTQ {
    tag "Human depletion (split mapped vs depleted FASTQ)"
    label 'process_low'

    conda "bioconda::samtools=1.17"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'library://daanjansen94/metatropics/samtools:v1.17' :
        'daanjansen94/samtools:v1.17' }"

    input:
    tuple val(meta), path(input)
    val(interleave)

    output:
    tuple val(meta), path("*_{1,2}.fastq.gz")      , optional:true, emit: fastq
    tuple val(meta), path("*_interleaved.fastq.gz"), optional:true, emit: interleaved
    tuple val(meta), path("*_singleton.fastq.gz")  , optional:true, emit: singleton
    tuple val(meta), path("*_other.fastq.gz")      , optional:true, emit: other
    path  "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def threads = (task.cpus as int) > 1 ? ((task.cpus as int) - 1) : 1
    if (args?.toString() =~ /(^|\s)-f\s|(^|\s)-F\s/) {
        error "SAMTOOLS_hFASTQ internally uses -f/-F; do not pass -f/-F via ext.args"
    }
    """
    set -euo pipefail

    # Depletion: emit mapped reads + depleted reads.

    if ${meta.single_end}; then
        # Mapped reads (primary alignments only).
        samtools fastq ${args} --threads ${threads} -F 0x904 \\
            -0 ${prefix}_mapped_1.fastq.gz \\
            -s /dev/null \\
            ${input}

        # Depleted reads (unmapped; drop secondary/supplementary).
        samtools fastq ${args} --threads ${threads} -f 4 -F 0x900 \\
            -0 ${prefix}_unmapped_1.fastq.gz \\
            -s /dev/null \\
            ${input}

        mv ${prefix}_mapped_1.fastq.gz ${prefix}_1.fastq.gz
        mv ${prefix}_unmapped_1.fastq.gz ${prefix}_other.fastq.gz
        # Drop empty gz placeholders.
        for f in ${prefix}_1.fastq.gz ${prefix}_other.fastq.gz; do
            if [ -f "\$f" ]; then
                perl -e 'exit((stat(shift))[7] <= 28 ? 0 : 1)' "\$f" && rm -f "\$f" || true
            fi
        done
    else
        if ${interleave}; then
            echo "ERROR: interleaved output is not supported for depletion mode" >&2
            exit 1
        fi

        samtools fastq ${args} --threads ${threads} -F 0x904 \\
            -1 ${prefix}_mapped_1.fastq.gz \\
            -2 ${prefix}_mapped_2.fastq.gz \\
            -s /dev/null \\
            -0 /dev/null \\
            ${input}

        samtools fastq ${args} --threads ${threads} -f 4 -F 0x900 \\
            -1 ${prefix}_unmapped_1.fastq.gz \\
            -2 ${prefix}_unmapped_2.fastq.gz \\
            -s /dev/null \\
            -0 /dev/null \\
            ${input}

        mv ${prefix}_mapped_1.fastq.gz ${prefix}_1.fastq.gz
        mv ${prefix}_mapped_2.fastq.gz ${prefix}_2.fastq.gz
        # Paired-end: publish depleted R1 only.
        mv ${prefix}_unmapped_1.fastq.gz ${prefix}_other.fastq.gz
        rm -f ${prefix}_unmapped_2.fastq.gz || true
        for f in ${prefix}_1.fastq.gz ${prefix}_2.fastq.gz ${prefix}_other.fastq.gz; do
            if [ -f "\$f" ]; then
                perl -e 'exit((stat(shift))[7] <= 28 ? 0 : 1)' "\$f" && rm -f "\$f" || true
            fi
        done
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
