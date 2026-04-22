process MEDAKA {
    tag "Variant calling"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'library://daanjansen94/metatropics/medaka:v1.4.3' :
        'daanjansen94/medaka:v1.4.3' }"

    input:
    tuple val(meta), path(reads), path(assembly)

    output:
    tuple val(meta), path("*.vcf"), optional: true, emit: assembly
    tuple val(meta), path("*.bam"), path("*.bai"), optional: true, emit: bamfiles
    path "versions.yml", emit: versions
    path "*.coverage.txt", optional: true, emit: coveragefiles

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}.${meta.virus}"
    """
    if [ -s $assembly ] && [ \$(grep -c ">" $assembly) -gt 0 ]; then
        medaka_haploid_variant \\
            -t $task.cpus \\
            $args \\
            -i $reads \\
            -r $assembly \\
            -o ./

        if [ \$? -eq 0 ]; then
            mv medaka.annotated.vcf ${prefix}.vcf
            mv calls_to_ref.bam ${prefix}.sorted.bam
            mv calls_to_ref.bam.bai ${prefix}.sorted.bam.bai
            rm -f medaka.vcf

            samtools depth -a ${prefix}.sorted.bam > ${prefix}.sorted.bam.coverage.txt
        else
            echo "Medaka processing failed for ${prefix}. Creating empty output files." >&2
            touch ${prefix}.vcf
            touch ${prefix}.sorted.bam
            touch ${prefix}.sorted.bam.bai
        fi
    else
        echo "Assembly file is empty or contains no sequences for ${prefix}. Creating empty output files." >&2
        touch ${prefix}.vcf
        touch ${prefix}.sorted.bam
        touch ${prefix}.sorted.bam.bai
    fi

    # Remove empty files
    find . -type f -empty -delete

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        medaka: \$(medaka --version 2>&1 | sed 's/medaka //g')
        samtools: \$(samtools --version | sed '1!d; s/samtools //')
    END_VERSIONS
    """
}
