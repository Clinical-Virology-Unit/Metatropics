process MEDAKA_VARIANTS {
    tag "Medaka variant calling"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'library://jansendaan94_v2/metatropics/medaka:2.0.0' :
        'daanjansen94/medaka:2.0.0' }"

    input:
    tuple val(meta), path(reads), path(assembly)

    output:
    tuple val(meta), path("*.vcf"), optional: true, emit: assembly
    tuple val(meta), path("*.medaka.filtered.vcf"), emit: filtered
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}.${meta.virus_slug}"
    def user_model = params.medaka_model ?: ''
    def hac_fallback = 'r1041_e82_400bps_hac_variant_v5.0.0'
    def min_qual = params.quality
    def min_dpsp = params.depth
    """
    write_skeleton_filtered_vcf() {
        # Medaka did not emit a usable annotated VCF. Downstream postprocessing still expects a valid VCF header.
        : > "${prefix}.medaka.filtered.vcf"
        echo '##fileformat=VCFv4.2' >> "${prefix}.medaka.filtered.vcf"
        if [ ! -s "${assembly}.fai" ]; then
            samtools faidx $assembly
        fi
        cut -f1,2 "${assembly}.fai" | while read -r chrom clen; do
            echo "##contig=<ID=\${chrom},length=\${clen}>" >> "${prefix}.medaka.filtered.vcf"
        done
        echo '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">' >> "${prefix}.medaka.filtered.vcf"
        echo '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE' >> "${prefix}.medaka.filtered.vcf"
    }

    if [ -s $assembly ] && [ \$(grep -c ">" $assembly) -gt 0 ]; then
        # Force CPU execution.
        # Medaka's device selection is based on torch CUDA detection; in some container/runtime
        # configurations it can think a GPU exists but CUDA isn't actually usable, causing
        # model loading to fail. For Metatropics we prefer a robust CPU-only run.
        export CUDA_VISIBLE_DEVICES=""

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
            # Filter Medaka's annotated VCF: drop variants with low confidence (QUAL) or low
            # spoa-realigned spanning depth (DPSP). SNP/indel split lives in postprocessing.
            bcftools view -e 'QUAL<${min_qual} || INFO/DPSP<${min_dpsp}' \\
                medaka.annotated.vcf > medaka.filtered.vcf
            cp medaka.filtered.vcf ${prefix}.medaka.filtered.vcf
        else
            echo "Medaka produced no annotated VCF for ${prefix}." >&2
            write_skeleton_filtered_vcf
        fi
    else
        echo "Assembly file is empty or contains no sequences for ${prefix}. Skipping Medaka." >&2
        write_skeleton_filtered_vcf
    fi

    # Drop Medaka intermediate VCFs; keep `${prefix}.medaka.filtered.vcf` for downstream.
    rm -f medaka.vcf medaka.sorted.vcf medaka.annotated.vcf medaka.filtered.vcf

    # Remove empty placeholder VCFs so the optional emit stays empty on failure.
    # Keep `${prefix}.medaka.filtered.vcf` even if it only contains a header (no variant rows).
    find . -maxdepth 1 -type f -name '*.vcf' ! -name '*.medaka.filtered.vcf' -empty -delete

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        medaka: \$(medaka --version 2>/dev/null | sed 's/medaka //')
        bcftools: \$(bcftools --version 2>/dev/null | sed -n '1s/^bcftools //p')
    END_VERSIONS
    """
}
