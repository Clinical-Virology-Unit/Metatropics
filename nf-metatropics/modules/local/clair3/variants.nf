process CLAIR3_VARIANTS {
    tag "Clair3 variant calling"
    label 'process_medium'

    // Use the Docker image for Singularity too (pulled via docker://) because the Sylabs library tag
    // may not exist for all architectures.
    container "${ workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer' ?
        'docker://daanjansen94/clair3:v2.0.1' :
        'daanjansen94/clair3:v2.0.1' }"

    input:
    tuple val(meta), path(bam), path(bai), path(ref_fasta), path(raw_fastq_gz)

    output:
    tuple val(meta), path("${meta.id}.${meta.virus_slug}.clair3.vcf.gz")    , emit: vcf_gz
    tuple val(meta), path("${meta.id}.${meta.virus_slug}.clair3.vcf.gz.tbi"), emit: vcf_tbi
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}.${meta.virus_slug}"
    def model_override = params.clair3_model ?: ''
    def min_mq = params.clair3_min_mq
    // Use a permissive coverage threshold for the caller itself.
    // Clair3 requires --min_coverage >= 2.
    // Depth cutoffs are enforced later during uniform post-processing and consensus.
    def min_cov = 2
    def min_qual = params.quality
    def threads = task.cpus ?: 8
    """
    set -euo pipefail

    # Stage inputs locally (Nextflow often uses symlinks).
    # Convert contig names containing `|` (and sometimes `:`) to `_` for Clair3/htslib,
    # in BOTH FASTA and BAM header (must match).
    cp -f ${ref_fasta} ref.fasta
    cp -f ${bam} in.bam
    cp -f ${bai} in.bam.bai

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
    CTG_NAME=\$(cut -f1 ref.pipefix.fasta.fai | head -n 1)

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

    # Autodetect Clair3 model from BAM header (Dorado/Guppy model tags) when possible.
    # If the BAM header doesn't contain basecaller/model tags, try the raw/fixed FASTQ used in the pipeline
    # (these headers typically retain the RG:Z:...r10.4.1_e8.2_400bps_hac@v5.0.0 tag). Prefer HAC.
    DETECTED_HAC=\$(samtools view -H in.bam | tr '\\t' '\\n' | grep -Eo 'r[0-9]+[^ ]*_hac_v[0-9]+' | head -n 1 || true)
    DETECTED_SUP=\$(samtools view -H in.bam | tr '\\t' '\\n' | grep -Eo 'r[0-9]+[^ ]*_sup_v[0-9]+' | head -n 1 || true)
    if [ -z "\${DETECTED_HAC}" ] && [ -z "\${DETECTED_SUP}" ] && [ -s "${raw_fastq_gz}" ]; then
        RAW_TAG=\$(zcat ${raw_fastq_gz} 2>/dev/null | head -n 40000 | grep -Eo 'r10\\.[0-9]+\\.[0-9]+_e[0-9]+\\.[0-9]+_400bps_(hac|sup)@v[0-9]+\\.[0-9]+\\.[0-9]+' | head -n 1 || true)
        if [ -n "\${RAW_TAG}" ]; then
            # Example in raw reads: r10.4.1_e8.2_400bps_hac@v5.0.0  →  r1041_e82_400bps_hac_v500
            NORMALIZED=\$(echo "\${RAW_TAG}" | sed -E 's/\\.//g; s/@v/_v/')
            DETECTED_HAC=\$(echo "\${NORMALIZED}" | grep -Eo 'r[0-9]+[^ ]*_hac_v[0-9]+' | head -n 1 || true)
            DETECTED_SUP=\$(echo "\${NORMALIZED}" | grep -Eo 'r[0-9]+[^ ]*_sup_v[0-9]+' | head -n 1 || true)
        fi
    fi
    if [ -n "${model_override}" ]; then
        MODEL_NAME="${model_override}"
        echo "[clair3] Using user-specified model override: \${MODEL_NAME}" >&2
    elif [ -n "\${DETECTED_HAC}" ] && [ -d "/opt/models/\${DETECTED_HAC}" ]; then
        MODEL_NAME="\${DETECTED_HAC}"
        echo "[clair3] Auto-detected HAC model from BAM header: \${MODEL_NAME}" >&2
    elif [ -n "\${DETECTED_SUP}" ] && [ -d "/opt/models/\${DETECTED_SUP}" ]; then
        MODEL_NAME="\${DETECTED_SUP}"
        echo "[clair3] Auto-detected SUP model from BAM header: \${MODEL_NAME}" >&2
    else
        # Default fallback: HAC (per your preference). Avoid *_with_mv unless you know your image supports it.
        MODEL_NAME="r1041_e82_400bps_hac_v500"
        echo "[clair3] Could not reliably detect model; using fallback: \${MODEL_NAME}" >&2
    fi

    /opt/bin/run_clair3.sh \\
      --bam_fn bam.pipefix.bam \\
      --ref_fn ref.pipefix.fasta \\
      --model_path "/opt/models/\${MODEL_NAME}" \\
      --threads ${threads} \\
      --platform ont \\
      --ctg_name "\${CTG_NAME}" \\
      --sample_name "${meta.id}" \\
      --qual ${min_qual} \\
      --min_mq ${min_mq} \\
      --min_coverage ${min_cov} \\
      --haploid_sensitive 1 \\
      --output ./

    # Clair3 outputs: merge_output.vcf.gz (+tbi) is the final merged callset.
    if [ -s merge_output.vcf.gz ] && [ ! -s merge_output.vcf.gz.tbi ]; then
        if command -v tabix >/dev/null 2>&1; then
            tabix -f -p vcf merge_output.vcf.gz
        elif command -v bcftools >/dev/null 2>&1; then
            bcftools index -f -t merge_output.vcf.gz
        else
            echo "[clair3] ERROR: merge_output.vcf.gz.tbi missing and neither tabix nor bcftools is available to create it." >&2
            exit 1
        fi
    fi
    cp -f merge_output.vcf.gz     ${prefix}.clair3.vcf.gz
    cp -f merge_output.vcf.gz.tbi ${prefix}.clair3.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clair3: \$(python3 /opt/bin/clair3.py --version 2>/dev/null || true)
        samtools: \$(samtools --version 2>/dev/null | head -n1 | sed 's/samtools //')
    END_VERSIONS
    """
}

