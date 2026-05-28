process CLAIR3_POSTPROCESSING {
    tag "Clair3 postprocessing"
    label 'process_low'

    // Clair3 image already includes python+pysam+samtools, so we can run postprocessing there.
    // Use the Docker image for Singularity too (pulled via docker://) because the Sylabs library tag
    // may not exist for all architectures.
    container "${ workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer' ?
        'docker://daanjansen94/clair3:v2.0.1' :
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

    # Clair3 VCF was generated using a contig-name fix (| and : -> _).
    # Apply the same transformation to the BAM and reference FASTA so allele recount and context lookup match.
    cp -f $ref_fasta ref.fasta
    cp -f $bam in.bam
    cp -f $bai in.bam.bai

    awk '
      BEGIN{OFS=""}
      /^>/{
        n=split(substr(\$0,2), a, " ")
        gsub(/[\\|:]/, "_", a[1])
        printf ">%s", a[1]
        for(i=2;i<=n;i++) printf " %s", a[i]
        printf "\\n"
        next
      }
      {print}
    ' ref.fasta > ref.pipefix.fasta
    samtools faidx ref.pipefix.fasta

    samtools view -H in.bam | awk '
      BEGIN{OFS="\\t"}
      /^@SQ/{
        for(i=1;i<=NF;i++){
          if(\$i ~ /^SN:/){
            sub(/^SN:/,"",\$i); gsub(/[\\|:]/,"_",\$i); \$i="SN:"\$i
          }
        }
      }
      {print}
    ' > bam.header.pipefix.sam
    samtools reheader bam.header.pipefix.sam in.bam > bam.pipefix.bam
    samtools index bam.pipefix.bam

    python ${projectDir}/bin/clair3.py \\
        --sample '${meta.id}' \\
        --virus '${meta.virus}' \\
        --bam bam.pipefix.bam \\
        --ref-fasta ref.pipefix.fasta \\
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

