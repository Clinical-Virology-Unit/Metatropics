/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Resolve Host keywords into concrete FASTA paths (do not mutate params)
def resolvedHosts = HostReferences.resolve(params, log, workflow.projectDir.parent)
def resolvedHumanHostFasta = params.Human_host_fasta ?: resolvedHosts.human ?: params.fasta
def resolvedOtherHostFasta = params.Other_host_fasta ?: resolvedHosts.other ?: params.host_fasta

// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, resolvedHumanHostFasta, resolvedOtherHostFasta ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK_METATROPICS } from './subworkflows/local/input_check_metatropics'
include { FIX } from './subworkflows/local/subfix_names'
include { HUMAN_MAPPING } from './subworkflows/local/human_mapping'
include { HOST_MAPPING } from './subworkflows/local/host_mapping'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { DORADO_ONT } from '../modules/local/dorado/ont'
include { DORADO_DEMULTIPLEXING } from '../modules/local/dorado/demultiplexing'
include { RAREFACTION		          } from '../modules/local/rarefaction/rarefaction'
include { FASTPLONG                   } from '../modules/local/fastplong/main'
include { NANOPLOT                    } from '../modules/nf-core/nanoplot/main'
include { VIRASIGN_CLASSIFICATION      } from '../modules/local/virasign/classification'
include { VIRASIGN_DB                  } from '../modules/local/virasign/prepare_db'
include { VIRASIGN_SUMMARY             } from '../modules/local/virasign/build_html'
include { METAMAPS_MAP                } from '../modules/local/metamaps/map'
include { METAMAPS_CLASSIFY           } from '../modules/local/metamaps/classify'
include { R_METAPLOT                  } from '../modules/local/r/metaplot'
include { REF_FASTA                   } from '../modules/local/ref_fasta'
include { SEQTK_SUBSEQ                } from '../modules/nf-core/seqtk/subseq/main'
include { REFFIX_FASTA                } from '../modules/local/reffix_fasta'
include { MEDAKA                      } from '../modules/nf-core/medaka/main'
include { ReadCount                   } from '../modules/local/reads/reads'
include { RCOVERAGE                   } from '../modules/local/rcoverage/rcoverage'
include { SAMTOOLS_COVERAGE           } from '../modules/nf-core/samtools/coverage/main'
include { IVAR_CONSENSUS              } from '../modules/nf-core/ivar/consensus/main'
include { HOMOPOLISH_POLISHING        } from '../modules/local/homopolish/polishing'
include { ADDING_DEPTH                } from '../modules/local/adding_depth'
include { FINAL_REPORT                } from '../modules/local/final_report'
include { CLEANUP		              } from '../modules/local/cleanup/cleanup'
include { CLEANUP_INTERMEDIATE		  } from '../modules/local/cleanup/cleanup_intermediate'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


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

    FASTPLONG(
        ch_reads_for_fastp
    )

    NANOPLOT(
         ch_fixed_reads
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

        // Barrier = count(per-sample JSON (final if present, otherwise unfiltered)).
        // VIRASIGN_SUMMARY reads `${params.outdir}/Classification/virasign` inside the container.
        VIRASIGN_SUMMARY(
            VIRASIGN_CLASSIFICATION.out.final_json
                .mix(VIRASIGN_CLASSIFICATION.out.unfiltered_json)
                .count()
        )
    }

    METAMAPS_MAP(
        readsForViralClassification
    )

    meta_with_othermeta = METAMAPS_MAP.out.metaclass.join(METAMAPS_MAP.out.otherclassmeta)
    meta_with_othermeta_with_metalength = meta_with_othermeta.join(METAMAPS_MAP.out.metalength)
    meta_with_othermeta_with_metalength_with_parameter = meta_with_othermeta_with_metalength.join(METAMAPS_MAP.out.metaparameters)

    METAMAPS_CLASSIFY(
        meta_with_othermeta_with_metalength_with_parameter
    )

    // Collect all outputs from METAMAPS_CLASSIFY
    ch_all_metamaps_classify = METAMAPS_CLASSIFY.out.classem.collect()
    .mix(METAMAPS_CLASSIFY.out.classem_original.collect())
    .mix(METAMAPS_CLASSIFY.out.classlength.collect())
    .mix(METAMAPS_CLASSIFY.out.classWIMP.collect())
    .mix(METAMAPS_CLASSIFY.out.classcov.collect())
    .collect()

    // Create a cleanup channel with a dummy file path
    ch_cleanup_done = Channel.value(file("${workDir}/.cleanup_done"))

    // Run intermediate CLEANUP only if Docker cleanup is enabled
    if (params.enable_docker_cleanup) {
    CLEANUP_INTERMEDIATE(ch_all_metamaps_classify)
    ch_cleanup_done = CLEANUP_INTERMEDIATE.out.cleanup_done
    } else {
    // Create a dummy file in the work directory
    file("${workDir}/.cleanup_done").text = "Cleanup not enabled"
    }

    // Continue with the rest of your workflow, using ch_cleanup_done to ensure cleanup has finished
    rmetaplot_ch = ((METAMAPS_MAP.out.metaclass.join(METAMAPS_CLASSIFY.out.classlength))
            .join(METAMAPS_CLASSIFY.out.classcov))
            .join(NANOPLOT.out.totalreads)
            .combine(ch_cleanup_done)

    R_METAPLOT(
    rmetaplot_ch.map { meta, classification_results, length_and_identities, contig_coverage, total_reads, cleanup_done ->
        tuple(meta, classification_results, length_and_identities, contig_coverage, total_reads, cleanup_done)
    }
    )

    reffasta_ch=(R_METAPLOT.out.reporttsv.join(METAMAPS_CLASSIFY.out.classem)).join(readsForViralClassification)

    REF_FASTA(
        reffasta_ch
    )

	fixingheader_ch = REF_FASTA.out.headereads.map { entry ->
    def meta = entry[0]
    def files = entry[1]
    
    if (files instanceof Path) {
        [[id: meta.id, single_end: meta.single_end], [files]]  // Single file case
    } else {
        [[id: meta.id, single_end: meta.single_end], files]    // Multiple files case
	}
	}

	fixiseqref_ch = REF_FASTA.out.seqref.map { entry ->
    def meta = entry[0]
    def files = entry[1]
    
    if (files instanceof Path) {
        [[id: meta.id, single_end: meta.single_end], [files]]  // Single file case
    } else {
        [[id: meta.id, single_end: meta.single_end], files]    // Multiple files case
	}
	}

	fixingallreads_ch = REF_FASTA.out.allreads.map { entry ->
    def meta = entry[0]
    def files = entry[1]
    
    if (files instanceof Path) {
        [[id: meta.id, single_end: meta.single_end], [files]]  // Single file case
    } else {
        [[id: meta.id, single_end: meta.single_end], files]    // Multiple files case
	}
	}

	// FlatMap function for headers
	headers_ch = fixingheader_ch.flatMap { entry ->
        def id = entry[0].id
        def singleEnd = entry[0].single_end
        entry[1].collect { virus ->
            [[id: id, single_end: singleEnd, virus: virus.getBaseName().replaceFirst(/.+\./,"")], "${virus}"]
        }
    }

	// FlatMap function for ref
	fasta_ch = fixiseqref_ch.flatMap { entry ->
        def id = entry[0].id
        def singleEnd = entry[0].single_end
        def tm = entry[1].size()
              entry[1].collect { virus ->
                [[id: id, single_end: singleEnd, virus: ((virus.getBaseName()).replaceFirst(/\.REF+/,"")).replaceFirst(/.+\./,"")],  "${virus}"]
            }
    }

	// FlatMap function for fastq
	fastq_ch = fixingallreads_ch.flatMap { entry ->
        def id = entry[0].id
        def singleEnd = entry[0].single_end
        entry[1].collect { virus ->
            [[id: id, single_end: singleEnd, virus: virus.getBaseName().replaceFirst(/.+\./,"")], "${virus}"]
        }
    }
    
	//Ending of the fix channels per pathogen.


    REFFIX_FASTA(
        fasta_ch
    )

    SEQTK_SUBSEQ(
        fastq_ch.join(headers_ch)
    )

    MEDAKA(
        SEQTK_SUBSEQ.out.sequences.join(REFFIX_FASTA.out.fixedseqref)
    )

    // Define the host_genome_status
    def host_genome_status = 'not_used'
    if (resolvedHumanHostFasta && resolvedOtherHostFasta) {
        host_genome_status = 'both'
    } else if (resolvedHumanHostFasta) {
        host_genome_status = 'human_only'
    } else if (resolvedOtherHostFasta) {
        host_genome_status = 'other_only'
    }

    // Call the ReadCount process
    ReadCount(
    params.outdir,
    MEDAKA.out.coveragefiles.collect(),
    host_genome_status
    )

   // Conditional RCOVERAGE process
   if (params.rcoverage_figure) {
    RCOVERAGE(
        MEDAKA.out.coveragefiles.collect()
    )
    ch_rcoverage_done = RCOVERAGE.out.collect() // Create a channel that signals RCOVERAGE is done
    } else {
    ch_rcoverage_done = Channel.empty() // Create an empty channel if RCOVERAGE is not run
    }

    SAMTOOLS_COVERAGE(
        MEDAKA.out.bamfiles
    )

    savempileup = false
    IVAR_CONSENSUS(
        MEDAKA.out.bamfiles.join(REFFIX_FASTA.out.fixedseqref),
        savempileup
    )

    HOMOPOLISH_POLISHING(
        IVAR_CONSENSUS.out.fasta.join(REFFIX_FASTA.out.fixedseqref)
    )

    group_virus_and_ref_ch = (HOMOPOLISH_POLISHING.out.polishconsensus).map { entry ->
        def id = entry[0].id
        def singleEnd = entry[0].single_end
        def virus = entry[0].virus
        //def fasta = entry[1],entry[2]
        [[virus: virus], entry[1]]
    }.groupTuple()//.view()

    covconstat_ch = (SAMTOOLS_COVERAGE.out.coverage
        .join(HOMOPOLISH_POLISHING.out.polishconsensus)
        .join(SAMTOOLS_COVERAGE.out.bamstats)
    ).map { entry ->
        // entry = [meta, coverage_txt, consensus, bamstats]
        [[id: entry[0].id, single_end: entry[0].single_end], entry[1], entry[2], entry[3]]
    }

    addingdepthin_ch = (covconstat_ch.combine(R_METAPLOT.out.reporttsv, by: 0)).map { entry ->
        // entry = [meta_simple, coverage_txt, consensus, bamstats, report_tsv]
        def id = entry[0].id
        def singleEnd = entry[0].single_end
        def virus = entry[1].getBaseName().replaceFirst(/.+\./,"")
        // Order matches ADDING_DEPTH input: meta, depth, consensus, bamstats, report
        [[id: id, single_end: singleEnd, virus: virus], entry[1], entry[2], entry[3], entry[4]]
    }

    ADDING_DEPTH(
        addingdepthin_ch
    )

    FINAL_REPORT(
        (ADDING_DEPTH.out.repdepth.map{it[1]}).collect()
    )


    ch_versions = ch_versions.mix(FASTPLONG.out.versions.first())
    ch_versions = ch_versions.mix(NANOPLOT.out.versions.first())
    if (params.run_virasign) {
        ch_versions = ch_versions.mix(VIRASIGN_DB.out.versions)
    }
    ch_versions = ch_versions.mix(METAMAPS_MAP.out.versions.first())
    ch_versions = ch_versions.mix(METAMAPS_CLASSIFY.out.versions.first())
    ch_versions = ch_versions.mix(R_METAPLOT.out.versions.first())
    ch_versions = ch_versions.mix(SEQTK_SUBSEQ.out.versions.first())
    ch_versions = ch_versions.mix(MEDAKA.out.versions.first())
    ch_versions = ch_versions.mix(SAMTOOLS_COVERAGE.out.versions.first())
    ch_versions = ch_versions.mix(IVAR_CONSENSUS.out.versions.first())
    ch_versions = ch_versions.mix(HOMOPOLISH_POLISHING.out.versions.first())
    if (resolvedHumanHostFasta) {
        ch_versions = ch_versions.mix(HUMAN_MAPPING.out.versionsmini)
        ch_versions = ch_versions.mix(HUMAN_MAPPING.out.versionssamsort)
        ch_versions = ch_versions.mix(HUMAN_MAPPING.out.versionssamfastq)
    }

    // Wait for RCOVERAGE to complete before running CUSTOM_DUMPSOFTWAREVERSIONS
    CUSTOM_DUMPSOFTWAREVERSIONS(
    ch_versions.unique().collectFile(name: 'collated_versions.yml'),
    ch_rcoverage_done // Add this channel as an input
    )

    // Run CLEANUP only if Docker cleanup is enabled
    if (params.enable_docker_cleanup) {
        CLEANUP(
            CUSTOM_DUMPSOFTWAREVERSIONS.out.versions,
            FINAL_REPORT.out.finalReport,
            ReadCount.out.read_counts_csv
        )
    }

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