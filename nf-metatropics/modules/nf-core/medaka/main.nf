process MEDAKA {
    tag "Variant calling"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'library://jansendaan94_v2/metatropics/medaka:2.0.0' :
        'daanjansen94/medaka:2.0.0' }"

    input:
    tuple val(meta), path(reads), path(assembly)

    output:
    tuple val(meta), path("*.vcf"), optional: true, emit: assembly
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}.${meta.virus}"
    def user_model = params.medaka_model ?: ''
    def hac_fallback = 'r1041_e82_400bps_hac_variant_v5.0.0'
    def min_qual = params.quality
    def min_dpsp = params.depth
    """
    if [ -s $assembly ] && [ \$(grep -c ">" $assembly) -gt 0 ]; then
        if [ -n "${user_model}" ]; then
            MEDAKA_MODEL="${user_model}"
            echo "Using user-specified Medaka model: \${MEDAKA_MODEL}" >&2
        else
            MEDAKA_MODEL=\$(medaka tools resolve_model --auto_model variant $reads 2>/dev/null || true)
            if [ -z "\${MEDAKA_MODEL}" ]; then
                MEDAKA_MODEL="${hac_fallback}"
                echo "No Dorado basecaller tag detected on $reads. Falling back to HAC model: \${MEDAKA_MODEL}" >&2
            else
                echo "Auto-detected Medaka model from Dorado tag: \${MEDAKA_MODEL}" >&2
            fi
        fi

        medaka_variant \\
            -t $task.cpus \\
            -m \${MEDAKA_MODEL} \\
            $args \\
            -i $reads \\
            -r $assembly \\
            -o ./

        if [ -s medaka.annotated.vcf ]; then
            # Filter Medaka's annotated VCF: drop variants with low confidence
            # (QUAL) or low spoa-realigned spanning depth (DPSP). Then split the
            # surviving variants into separate SNP-only and indel-only files for
            # easier downstream inspection.
            bcftools view -e 'QUAL<${min_qual} || INFO/DPSP<${min_dpsp}' \\
                medaka.annotated.vcf > medaka.filtered.vcf
            bcftools view -v snps   medaka.filtered.vcf > ${prefix}_snps.vcf
            bcftools view -v indels medaka.filtered.vcf > ${prefix}_indel.vcf
        else
            echo "Medaka produced no annotated VCF for ${prefix}." >&2
            touch ${prefix}_snps.vcf
            touch ${prefix}_indel.vcf
        fi
    else
        echo "Assembly file is empty or contains no sequences for ${prefix}. Skipping Medaka." >&2
        touch ${prefix}_snps.vcf
        touch ${prefix}_indel.vcf
    fi

    # Drop Medaka's own intermediate VCFs so only the per-virus _snps/_indel
    # files are emitted.
    rm -f medaka.vcf medaka.sorted.vcf medaka.annotated.vcf medaka.filtered.vcf

    # Remove empty placeholder VCFs so the optional emit stays empty on failure.
    find . -maxdepth 1 -name '*.vcf' -type f -empty -delete

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        medaka: \$(medaka --version 2>/dev/null | sed 's/medaka //')
        bcftools: \$(bcftools --version 2>/dev/null | sed -n '1s/^bcftools //p')
    END_VERSIONS
    """
}
