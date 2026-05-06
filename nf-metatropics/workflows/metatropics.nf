/* Validate inputs */

import groovy.json.JsonSlurper

// Resolve Host keywords into concrete FASTA paths (do not mutate params)
def resolvedHosts = HostReferences.resolve(params, log, workflow.projectDir.parent)
def resolvedHumanHostFasta = params.Human_host_fasta ?: resolvedHosts.human ?: params.fasta
def resolvedOtherHostFasta = params.Other_host_fasta ?: resolvedHosts.other ?: params.host_fasta

// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, resolvedHumanHostFasta, resolvedOtherHostFasta ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

/* Imports: local subworkflows */
include { INPUT_CHECK_METATROPICS } from './subworkflows/local/input_check_metatropics'
include { FIX } from './subworkflows/local/subfix_names'
include { HUMAN_MAPPING } from './subworkflows/local/human_mapping'
include { HOST_MAPPING } from './subworkflows/local/host_mapping'

/* Imports: modules */
include { CUSTOM_DUMPSOFTWAREVERSIONS as SOFTWARE_VERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { DORADO_ONT } from '../modules/local/dorado/ont'
include { DORADO_DEMULTIPLEXING } from '../modules/local/dorado/demultiplexing'
include { RAREFACTION		          } from '../modules/local/rarefaction/rarefaction'
include { FASTPLONG                   } from '../modules/local/fastplong/main'
include { NANOPLOT                    } from '../modules/nf-core/nanoplot/main'
include { VIRASIGN_CLASSIFICATION      } from '../modules/local/virasign/classification'
include { VIRASIGN_DB                  } from '../modules/local/virasign/prepare_db'
include { VIRASIGN_SUMMARY as METATROPICS_SUMMARY } from '../modules/local/virasign/build_html'
include { CLAIR3_VARIANTS              } from '../modules/local/clair3/variants'
include { CLAIR3_POSTPROCESSING        } from '../modules/local/clair3/postprocessing'
include { CONSENSUS_BCFTOOLS           } from '../modules/local/bcftools/consensus'
include { ReadCount                   } from '../modules/local/reads/reads'

/* Main workflow */


workflow METATROPICS {

    ch_versions = Channel.empty()
    def ch_fixed_reads

    // Auto-detect mode without mutating `params` (Nextflow may ignore re-assignments).
    def doBasecall = (params.basecall != null) ? (params.basecall as boolean) : (params.input_dir != null)

    // Validate input parameters (must run after -params-file is loaded)
    WorkflowMetatropics.initialise(params, log)

    INPUT_CHECK_METATROPICS{
        ch_input
        //ch_input2
    }

    if (doBasecall) {
        if (params.input_dir==null) { exit 1, 'POD5 input dir not specified!'}
        if (params.input==null) { exit 1, 'Sample sheet not specified!'}
        
        ch_sample = INPUT_CHECK_METATROPICS.out.reads.map{tuple(it[1],it[0])}
        ch_sample_sheet = Channel.fromPath(params.input, checkIfExists: true)

        inPOD5 = channel.fromPath(params.input_dir)

        DORADO_ONT(
            inPOD5
        )

        DORADO_DEMULTIPLEXING(
            DORADO_ONT.out.basecalling_ch
        )
        
        ch_barcode = DORADO_DEMULTIPLEXING.out.demultiplexed_fastq.flatten().map{file -> tuple(file.simpleName, file)}
        ch_sample_barcode = ch_sample.join(ch_barcode)

        FIX(
            ch_sample_barcode
        )
        ch_fixed_reads = FIX.out.reads

        ch_versions = ch_versions.mix(DORADO_ONT.out.versions)
        ch_versions = ch_versions.mix(DORADO_DEMULTIPLEXING.out.versions)
    }
    else {
        ch_sample = INPUT_CHECK_METATROPICS.out.reads.map{tuple(it[1].replaceFirst(/\/.+\//,""),it[0],it[1])}

        FIX(
            ch_sample
        )
        ch_fixed_reads = FIX.out.reads
    }

   // Conditional execution of RAREFACTION
   def ch_reads_for_fastp
   if (params.perform_rarefaction) {
    RAREFACTION(
        ch_fixed_reads,
        params.perform_rarefaction,
        params.target_bases
    )
        ch_reads_for_fastp = RAREFACTION.out.rarefied_reads
    } else {
        ch_reads_for_fastp = ch_fixed_reads
    }

    NANOPLOT(
         ch_fixed_reads
     )

    FASTPLONG(
        ch_reads_for_fastp
    )

    def readsAfterHuman = FASTPLONG.out.reads
    def readsForViralClassification

    if (resolvedHumanHostFasta) {
        HUMAN_MAPPING(
            FASTPLONG.out.reads,
            resolvedHumanHostFasta
        )
        readsAfterHuman = HUMAN_MAPPING.out.humanout
    }

    if (resolvedOtherHostFasta) {
        HOST_MAPPING(
            readsAfterHuman,
            resolvedOtherHostFasta
        )
        readsForViralClassification = HOST_MAPPING.out.hostout
    } else {
        readsForViralClassification = readsAfterHuman
    }

    // Depletion mode for ReadCount / readcount.py: not_used | human_only | other_only | both
    def host_genome_status = 'not_used'
    if (resolvedHumanHostFasta && resolvedOtherHostFasta) {
        host_genome_status = 'both'
    } else if (resolvedHumanHostFasta) {
        host_genome_status = 'human_only'
    } else if (resolvedOtherHostFasta) {
        host_genome_status = 'other_only'
    }

    // ── Virasign (phase 1): prepare DB once, then run per-sample with --no-html ──
    if (params.run_virasign) {
        // Shared on-host results tree, isolated per virasign_database to avoid mixing results
        // across parameter changes when running with `-resume`.
        def rawDbArg = params.virasign_database?.toString()?.trim()
        def effectiveDbArg = rawDbArg ?: 'RVDB'
        def virasignDbLabel = effectiveDbArg.replaceAll(/[^A-Za-z0-9._-]+/, '_')
        def virasignResultsRoot = file("${params.outdir}/Classification/virasign/${virasignDbLabel}")
        virasignResultsRoot.mkdirs()

        // Prepare DB once (prevents parallel workers from double-downloading).
        VIRASIGN_DB()

        // Per-sample Virasign (-o publish in work/), then merge into outdir; HTML pass reads that same tree.
        virasign_db_ready = VIRASIGN_DB.out.ready
        VIRASIGN_CLASSIFICATION(virasign_db_ready, readsForViralClassification)

        // Build per-(sample,accession) tuples directly from Virasign-produced JSONs.
        // IMPORTANT: Do not glob the outdir with `checkIfExists: true` here, because the files
        // only exist after VIRASIGN_CLASSIFICATION completes.
        ch_virasign_confident = VIRASIGN_CLASSIFICATION.out.final_json.flatMap { jsonFile ->
            def sampleDir = jsonFile.parent
            def sampleLabel = sampleDir.getBaseName()
            def sampleId = sampleLabel
                .replaceFirst(/\\.fastp$/, '')
                .replaceFirst(/(_T\\d+_other_T\\d+|_other_T\\d+|_T\\d+_other|_T\\d+|_other)$/, '')

            def hits = (List) (new JsonSlurper().parse(jsonFile))
            hits.collect { hit ->
                def acc = hit.accession?.toString()
                def rawSp = (hit.organism ?: hit.viral_species ?: '')?.toString()?.trim()
                if ((!rawSp || rawSp.isEmpty()) && hit.description) {
                    rawSp = hit.description.toString().trim().take(120)
                }
                def spSlug = (rawSp && !rawSp.isEmpty()) ? rawSp.replaceAll(/[^A-Za-z0-9._-]+/, '_').replaceAll(/^_+|_+$/, '').replaceAll(/_+/, '_') : ''
                def virusSlug = (spSlug && !spSlug.isEmpty()) ? "${acc}_${spSlug}" : acc
                def accDir = file("${sampleDir}/${acc}")
                def meta2 = [
                    id          : sampleId,
                    single_end : true,
                    virus       : acc,
                    virus_slug  : virusSlug,
                    species_slug: spSlug,
                ]
                [
                    meta2,
                    file("${accDir}/${acc}.bam"),
                    file("${accDir}/${acc}.bam.bai"),
                    file("${accDir}/${acc}.fasta"),
                    file("${accDir}/mread.fastq.gz")
                ]
            }
        }
    }

    def ch_readcount_barrier
    if (params.run_virasign) {
        ch_readcount_barrier = VIRASIGN_CLASSIFICATION.out.final_json
            .mix(VIRASIGN_CLASSIFICATION.out.unfiltered_json)
            .count()
    } else {
        ch_readcount_barrier = readsForViralClassification.count()
    }
    def ch_readcount_in = ch_readcount_barrier.map { b -> tuple(params.outdir, b, host_genome_status) }
    ReadCount( ch_readcount_in )


    if (!params.run_virasign) {
        exit 1, "This pipeline configuration requires 'run_virasign: true' to generate per-virus inputs for variant calling."
    }

    // Provide Clair3 with the fixed/raw reads used earlier in the pipeline so we can
    // autodetect the Dorado/Guppy model from the FASTQ header (when the BAM header lacks it).
    def ch_raw_reads_for_model = ch_fixed_reads.map { meta, fq -> [ meta.id, fq ] }

    def ch_clair3_in = ch_virasign_confident
        .map { meta, bam, bai, ref, virasign_reads -> [ meta.id, meta, bam, bai, ref ] }
        .join(ch_raw_reads_for_model, by: 0)
        .map { sample, meta, bam, bai, ref, raw_fq -> [ meta, bam, bai, ref, raw_fq ] }

    CLAIR3_VARIANTS( ch_clair3_in )

    def ch_clair3_uniform_in = ch_virasign_confident
        .map { meta, bam, bai, ref, reads -> [ meta, bam, bai, ref ] }
        .join(CLAIR3_VARIANTS.out.vcf_gz.join(CLAIR3_VARIANTS.out.vcf_tbi, by: 0), by: 0)
        .map { meta, bam, bai, ref, vcfgz, tbi -> [ meta, vcfgz, tbi, bam, bai, ref ] }

    CLAIR3_POSTPROCESSING( ch_clair3_uniform_in )

    def ch_consensus_in = CLAIR3_POSTPROCESSING.out.vcf
        .join(ch_virasign_confident.map { meta, bam, bai, ref, reads -> [ meta, ref ] }, by: 0)

    CONSENSUS_BCFTOOLS( ch_consensus_in )

    // Build final Metatropics summary only after draft consensuses are available.
    // (Consensus-derived breadth metrics are computed from these FASTAs.)
    METATROPICS_SUMMARY(
        CONSENSUS_BCFTOOLS.out.fasta.count()
    )


    ch_versions = ch_versions.mix(FASTPLONG.out.versions.first())
    ch_versions = ch_versions.mix(NANOPLOT.out.versions.first())
    if (params.run_virasign) {
        ch_versions = ch_versions.mix(VIRASIGN_DB.out.versions)
    }
    ch_versions = ch_versions.mix(CLAIR3_VARIANTS.out.versions.first())
    ch_versions = ch_versions.mix(CLAIR3_POSTPROCESSING.out.versions.first())
    ch_versions = ch_versions.mix(CONSENSUS_BCFTOOLS.out.versions.first())
    if (resolvedHumanHostFasta) {
        ch_versions = ch_versions.mix(HUMAN_MAPPING.out.versionsmini)
        ch_versions = ch_versions.mix(HUMAN_MAPPING.out.versionssamsort)
        ch_versions = ch_versions.mix(HUMAN_MAPPING.out.versionssamfastq)
    }

    SOFTWARE_VERSIONS(
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    NfcoreTemplate.summary(workflow, params, log)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/