process CONSENSUS_BCFTOOLS {
    tag "Consensus (bcftools, major tier)"
    label 'process_medium'

    // Use Docker image URI for Singularity/Apptainer pulls (Sylabs library tag may not exist for amd64).
    container "${ workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer' ?
        'docker://daanjansen94/bcftools:1.23.1' :
        'daanjansen94/bcftools:1.23.1' }"

    input:
    tuple val(meta), path(uniform_vcf), path(ref_fasta), path(depth_mask_bed)

    output:
    tuple val(meta), path("${meta.id}.${meta.virus_slug}.consensus.fasta"), emit: fasta
    tuple val(meta), path("${meta.id}.${meta.virus_slug}.consensus.metrics.json"), emit: metrics
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}.${meta.virus_slug}"
    def outFasta = "${meta.id}.${meta.virus_slug}.consensus.fasta"
    def outMetrics = "${meta.id}.${meta.virus_slug}.consensus.metrics.json"
    def consensus_vaf = (params.agreement != null) ? (params.agreement as Double) : 0.7d
    def min_depth = (params.depth != null) ? (params.depth as Integer) : 25
    def spp = (meta.species_slug ?: '').toString().trim()
    def fastaHdr = spp ? ">${meta.id}_${spp}" : ">${meta.id}"
    """
    set -euo pipefail

    # Virasign references/BAMs can contain `|` (and sometimes `:`) in contig names.
    # htslib/bcftools treats ':' as special (region syntax), so normalise BOTH to '_' everywhere.
    # Pipe-fix reference and mask so bcftools consensus sees the same contig names as the VCF.
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
    ' $ref_fasta > ref.pipefix.fasta

    awk 'BEGIN{OFS="\\t"} { gsub(/[\\|:]/, "_", \$1); print }' $depth_mask_bed > mask.pipefix.bed

    # Keep only "major" variants, but enforce a stricter VAF threshold for consensus assembly.
    # The tiered VCF (uniform_vcf) can keep a permissive major tier (e.g. 0.2),
    # while the consensus uses `params.agreement` (default 0.7).
    bcftools view -i "INFO/TIER=\\"major\\" && INFO/VAF>=${consensus_vaf}" -Oz -o ${prefix}.major.vcf.gz $uniform_vcf
    tabix -p vcf ${prefix}.major.vcf.gz

    # Apply externally-computed depth mask (depth < params.depth) to produce N-masked consensus.
    bcftools consensus -f ref.pipefix.fasta -m mask.pipefix.bed -o raw_consensus.fasta ${prefix}.major.vcf.gz
    awk 'NR==1{print "${fastaHdr}"; next} {print}' raw_consensus.fasta > ${outFasta}
    rm -f raw_consensus.fasta

    # Compute consensus breadth in reference coordinates so it:
    # - counts deletions as "determined" (they're still a call on reference positions)
    # - cannot exceed 100% (insertions do not increase reference length)
    REF_LEN=\$(awk '!/^>/{sum+=length(\$0)} END{printf(\"%d\",sum)}' ref.pipefix.fasta)
    MASKED=\$(awk '{sum+=\$3-\$2} END{printf(\"%d\",sum)}' mask.pipefix.bed)
    if [ -z "\$REF_LEN" ]; then REF_LEN=0; fi
    if [ -z "\$MASKED" ]; then MASKED=0; fi
    if [ "\$REF_LEN" -eq 0 ]; then BREADTH="0.0"; else BREADTH=\$(awk -v r="\$REF_LEN" -v m="\$MASKED" 'BEGIN{printf(\"%.4f\", (100.0*(r-m))/r)}'); fi
    printf '{\"ref_len\":%s,\"masked_bases_lt_min_depth\":%s,\"min_depth\":%s,\"consensus_breadth_pct\":%s}\\n' "\$REF_LEN" "\$MASKED" "${min_depth}" "\$BREADTH" > "${outMetrics}"

    rm -f ref.pipefix.fasta mask.pipefix.bed ${prefix}.major.vcf.gz ${prefix}.major.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>/dev/null | sed -n '1s/^bcftools //p')
        tabix: \$(tabix --version 2>/dev/null | head -n1 | sed -n 's/^tabix (htslib) //p')
    END_VERSIONS
    """
}

