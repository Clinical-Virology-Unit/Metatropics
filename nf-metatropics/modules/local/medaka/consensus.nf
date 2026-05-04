process MEDAKA_CONSENSUS_BCFTOOLS {
    tag "Consensus (bcftools, major tier)"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'library://jansendaan94_v2/metatropics/medaka:2.0.0' :
        'daanjansen94/medaka:2.0.0' }"

    input:
    tuple val(meta), path(uniform_vcf), path(ref_fasta)

    output:
    tuple val(meta), path("${meta.id}.${meta.virus_slug}.consensus.fasta"), emit: fasta
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}.${meta.virus_slug}"
    def outFasta = "${meta.id}.${meta.virus_slug}.consensus.fasta"
    def spp = (meta.species_slug ?: '').toString().trim()
    def fastaHdr = spp ? ">${meta.id}_${spp}" : ">${meta.id}"
    """
    set -euo pipefail

    samtools faidx $ref_fasta

    bcftools view -i 'INFO/TIER="major"' -Oz -o ${prefix}.major.vcf.gz $uniform_vcf
    tabix -p vcf ${prefix}.major.vcf.gz

    bcftools consensus -f $ref_fasta -o raw_consensus.fasta ${prefix}.major.vcf.gz
    awk 'NR==1{print "${fastaHdr}"; next} {print}' raw_consensus.fasta > ${outFasta}
    rm -f raw_consensus.fasta

    rm -f ${prefix}.major.vcf.gz ${prefix}.major.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>/dev/null | sed -n '1s/^bcftools //p')
        samtools: \$(samtools --version 2>/dev/null | head -n1 | sed 's/samtools //')
    END_VERSIONS
    """
}
