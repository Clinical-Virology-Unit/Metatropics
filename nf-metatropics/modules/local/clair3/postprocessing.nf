process CLAIR3_POSTPROCESSING {
    tag "Clair3 postprocessing"
    label 'process_low'

    // Clair3 image already includes python+pysam+samtools, so we can run postprocessing there.
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'library://jansendaan94_v2/metatropics/clair3:v2.0.1':
        'daanjansen94/clair3:v2.0.1' }"

    input:
    tuple val(meta), path(clair3_vcf_gz), path(clair3_vcf_tbi), path(bam), path(bai), path(ref_fasta)

    output:
    tuple val(meta), path("*.variants.filtered.vcf")   , emit: vcf
    tuple val(meta), path("*.variants.html")           , emit: html
    tuple val(meta), path("*.variants.unfiltered.vcf"), emit: variants_unfiltered
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}.${meta.virus_slug}"
    """
    set -euo pipefail

    python ${projectDir}/bin/clair3.py \\
        --sample '${meta.id}' \\
        --virus '${meta.virus}' \\
        --bam $bam \\
        --ref-fasta $ref_fasta \\
        --clair3-vcf $clair3_vcf_gz \\
        --out-prefix '${prefix}' \\
        --min-qual ${params.quality} \\
        --min-bq ${params.clair3_min_bq} \\
        --min-mq ${params.clair3_min_mq} \\
        --min-dp ${params.depth} \\
        --min-alt-reads ${params.clair3_min_alt_reads} \\
        --major-vaf ${params.major_vaf} \\
        --minor-vaf-min ${params.minor_vaf_min} \\
        --minor-vaf-max ${params.minor_vaf_max} \\
        --min-sb-pvalue ${params.min_sb_pvalue} \\
        --sb-min-alt-strand ${params.sb_min_alt_strand}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        bcftools: \$(bcftools --version 2>/dev/null | sed -n '1s/^bcftools //p')
        samtools: \$(samtools --version 2>/dev/null | head -n1 | sed 's/samtools //')
    END_VERSIONS
    """
}

