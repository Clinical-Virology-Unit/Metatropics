process MEDAKA_POSTPROCESSING {
    tag "Medaka postprocessing"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'library://jansendaan94_v2/metatropics/medaka:2.0.0' :
        'daanjansen94/medaka:2.0.0' }"

    input:
    tuple val(meta), path(medaka_filtered_vcf), path(bam), path(bai), path(ref_fasta)

    output:
    tuple val(meta), path("*.variants.uniform.vcf.gz")     , emit: vcf
    tuple val(meta), path("*.variants.uniform.vcf.gz.tbi"), emit: vcf_index
    tuple val(meta), path("*_snps.vcf")                   , emit: snps
    tuple val(meta), path("*_indel.vcf")                  , emit: indels
    tuple val(meta), path("*.variants.html")              , emit: html
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}.${meta.virus_slug}"
    """
    set -euo pipefail

    python ${projectDir}/bin/medaka.py \\
        --sample '${meta.id}' \\
        --virus '${meta.virus}' \\
        --bam $bam \\
        --ref-fasta $ref_fasta \\
        --medaka-filtered-vcf $medaka_filtered_vcf \\
        --out-prefix '${prefix}' \\
        --min-qual ${params.quality} \\
        --min-bq ${params.medaka_min_bq} \\
        --min-dp ${params.depth} \\
        --min-alt-reads ${params.medaka_min_alt_reads} \\
        --major-vaf ${params.medaka_major_vaf} \\
        --minor-vaf-min ${params.medaka_minor_vaf_min} \\
        --minor-vaf-max ${params.medaka_minor_vaf_max} \\
        --min-sb-pvalue ${params.medaka_min_sb_pvalue}

    bgzip -c ${prefix}.variants.uniform.vcf > ${prefix}.variants.uniform.vcf.gz
    tabix -p vcf ${prefix}.variants.uniform.vcf.gz
    rm -f ${prefix}.variants.uniform.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        medaka: \$(medaka --version 2>/dev/null | sed 's/medaka //')
        bcftools: \$(bcftools --version 2>/dev/null | sed -n '1s/^bcftools //p')
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """
}
