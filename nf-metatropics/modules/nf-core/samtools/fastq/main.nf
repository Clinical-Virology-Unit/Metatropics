process SAMTOOLS_hFASTQ {
    tag "Human depletion (keep depleted FASTQ only)"
    label 'process_low'

    conda "bioconda::samtools=1.17"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'library://daanjansen94/metatropics/samtools:v1.17' :
        'daanjansen94/samtools:v1.17' }"

    input:
    tuple val(meta), path(input)
    val(interleave)

    output:
    // Mapped human reads are intentionally NOT written (privacy).
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

    # Human depletion: emit ONLY unmapped (depleted) reads.
    # Do not extract mapped human reads to FASTQ (private / PHI risk).

    if ${meta.single_end}; then
        samtools fastq ${args} --threads ${threads} -f 4 -F 0x900 \\
            -0 ${prefix}_other.fastq.gz \\
            -s /dev/null \\
            ${input}

        if [ -f "${prefix}_other.fastq.gz" ]; then
            perl -e 'exit((stat(shift))[7] <= 28 ? 0 : 1)' "${prefix}_other.fastq.gz" && rm -f "${prefix}_other.fastq.gz" || true
        fi
    else
        if ${interleave}; then
            echo "ERROR: interleaved output is not supported for depletion mode" >&2
            exit 1
        fi

        samtools fastq ${args} --threads ${threads} -f 4 -F 0x900 \\
            -1 ${prefix}_other.fastq.gz \\
            -2 /dev/null \\
            -s /dev/null \\
            -0 /dev/null \\
            ${input}

        if [ -f "${prefix}_other.fastq.gz" ]; then
            perl -e 'exit((stat(shift))[7] <= 28 ? 0 : 1)' "${prefix}_other.fastq.gz" && rm -f "${prefix}_other.fastq.gz" || true
        fi
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
