process CONSENSUS_BCFTOOLS {
    tag "Consensus (bcftools, major tier)"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'library://daanjansen94/metatropics/bcftools:1.23.1' :
        'daanjansen94/bcftools:1.23.1' }"

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
    def consensus_vaf = (params.agreement != null) ? (params.agreement as Double) : 0.7d
    def spp = (meta.species_slug ?: '').toString().trim()
    def fastaHdr = spp ? ">${meta.id}_${spp}" : ">${meta.id}"
    """
    set -euo pipefail

    # Optional: index reference if samtools is available in the container.
    # (The bcftools biocontainer ships bcftools+tabix but not samtools.)
    if command -v samtools >/dev/null 2>&1; then
        samtools faidx $ref_fasta
    fi

    # Keep only "major" variants, but enforce a stricter VAF threshold for consensus assembly.
    # The tiered VCF (uniform_vcf) can keep a permissive major tier (e.g. 0.2),
    # while the consensus uses `params.agreement` (default 0.7).
    bcftools view -i "INFO/TIER=\\"major\\" && INFO/VAF>=${consensus_vaf}" -Oz -o ${prefix}.major.vcf.gz $uniform_vcf
    tabix -p vcf ${prefix}.major.vcf.gz

    bcftools consensus -f $ref_fasta -o raw_consensus.fasta ${prefix}.major.vcf.gz
    awk 'NR==1{print "${fastaHdr}"; next} {print}' raw_consensus.fasta > ${outFasta}
    rm -f raw_consensus.fasta

    rm -f ${prefix}.major.vcf.gz ${prefix}.major.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>/dev/null | sed -n '1s/^bcftools //p')
        tabix: \$(tabix --version 2>/dev/null | head -n1 | sed -n 's/^tabix (htslib) //p')
    END_VERSIONS
    """
}

