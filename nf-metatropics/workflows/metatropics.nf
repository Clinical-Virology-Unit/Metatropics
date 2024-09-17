
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowMetatropics.initialise(params, log)

// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config, params.fasta ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check'
include { INPUT_CHECK_METATROPICS } from '../subworkflows/local/input_check_metatropics'
include { FIX } from '../subworkflows/local/subfix_names'
include { HUMAN_MAPPING } from '../subworkflows/local/human_mapping'
include { HOST_MAPPING } from '../subworkflows/local/host_mapping'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
//include { FASTQC                      } from '../modules/nf-core/fastqc/main'
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { GUPPY_ONT                   } from '../modules/local/guppy/ont'
include { GUPPYDEMULTI_DEMULTIPLEXING } from '../modules/local/guppydemulti/demultiplexing'
include { RAREFACTION		      } from '../modules/local/rarefaction/rarefaction'
include { FASTP                       } from '../modules/nf-core/fastp/main'
include { NANOPLOT                    } from '../modules/nf-core/nanoplot/main'
include { METAMAPS_MAP                } from '../modules/local/metamaps/map'
include { METAMAPS_CLASSIFY           } from '../modules/local/metamaps/classify'
include { R_METAPLOT                  } from '../modules/local/r/metaplot'
include { KRONA_KRONADB               } from '../modules/nf-core/krona/kronadb/main'
include { KRONA_KTIMPORTTAXONOMY      } from '../modules/nf-core/krona/ktimporttaxonomy/main'
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
include { BAM_READCOUNT               } from '../modules/local/bam/readcount'
//include { MAFFT_ALIGN                 } from '../modules/local/mafft/align'
//include { SNIPIT_SNPPLOT              } from '../modules/local/snipit/snpplot'
//include { SNP_COMPARE                 } from '../modules/local/snp/compare'
//include { MAFFT_ALIGN as MAFFT_TWO    } from '../modules/local/mafft/align'
//include { SNIPIT_SNPPLOT as SNIPIT_TWO } from '../modules/local/snipit/snpplot'
include { CLEANUP		              } from '../modules/local/cleanup/cleanup'
include { CLEANUP_INTERMEDIATE		  } from '../modules/local/cleanup/cleanup_intermediate'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


// Info required for completion email and summary
def multiqc_report = []

workflow METATROPICS {

    ch_versions = Channel.empty()

    INPUT_CHECK_METATROPICS{
        ch_input
        //ch_input2
    }

    if(params.basecall==true){
        if (params.input_dir==null) { exit 1, 'FAST5 input dir not specified!'}
        ch_sample = INPUT_CHECK_METATROPICS.out.reads.map{tuple(it[1],it[0])}

        inFast5 = channel.fromPath(params.input_dir)

        GUPPY_ONT(
            inFast5
        )

        GUPPYDEMULTI_DEMULTIPLEXING(
            GUPPY_ONT.out.basecalling_ch
        )

        ch_barcode = GUPPYDEMULTI_DEMULTIPLEXING.out.barcodeReads.flatten().map{file -> tuple(file.simpleName, file)}
        ch_sample_barcode = ch_sample.join(ch_barcode)

        FIX(
            ch_sample_barcode
        )

        ch_versions = ch_versions.mix(GUPPY_ONT.out.versions)
        ch_versions = ch_versions.mix(GUPPYDEMULTI_DEMULTIPLEXING.out.versions)
    }
    else if(params.basecall==false){
        ch_sample = INPUT_CHECK_METATROPICS.out.reads.map{tuple(it[1].replaceFirst(/\/.+\//,""),it[0],it[1])}

        FIX(
            ch_sample
        )
    }

   // Define parameters for rarefaction
   params.perform_rarefaction = false
   params.target_bases = 1000000000 // Default value, can be overridden in the submission file

   // Conditional execution of RAREFACTION
   if (params.perform_rarefaction) {
    RAREFACTION(
        FIX.out.reads,
        params.perform_rarefaction,
        params.target_bases
    )
        ch_reads_for_fastp = RAREFACTION.out.rarefied_reads
    } else {
        ch_reads_for_fastp = FIX.out.reads
    }

    fastp_save_trimmed_fail = false
    FASTP(
        ch_reads_for_fastp,
        [],
        fastp_save_trimmed_fail,
        []
    )

    NANOPLOT(
         FIX.out.reads
     )

    HUMAN_MAPPING(
        FASTP.out.reads
    )

    if (params.host_fasta) {
        HOST_MAPPING(
            HUMAN_MAPPING.out.humanout
        )
        readsForMetamaps = HOST_MAPPING.out.hostout
    } else {
        readsForMetamaps = HUMAN_MAPPING.out.humanout
    }

    METAMAPS_MAP(
        readsForMetamaps
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
    .mix(METAMAPS_CLASSIFY.out.classkrona.collect())
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
    //KRONA_KRONADB();

    //KRONA_KTIMPORTTAXONOMY(
    //    METAMAPS_CLASSIFY.out.classkrona,
    //    KRONA_KRONADB.out.db
    //)

    reffasta_ch=(R_METAPLOT.out.reporttsv.join(METAMAPS_CLASSIFY.out.classem)).join(readsForMetamaps)

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
    def host_genome_status = params.host_fasta ? 'used' : 'not_used'

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

    covcon_ch = (SAMTOOLS_COVERAGE.out.coverage.join(HOMOPOLISH_POLISHING.out.polishconsensus)).map { entry ->
    [[id: entry[0].id, single_end: entry[0].single_end], entry[1], entry[2]]
    }

    addingdepthin_ch = (covcon_ch.combine(R_METAPLOT.out.reporttsv, by: 0)).map { entry ->
        def id = entry[0].id
        def singleEnd = entry[0].single_end
        def virus = entry[1].getBaseName().replaceFirst(/.+\./,"")
        [[id: id, single_end: singleEnd, virus: virus], entry[1], entry[2], entry[3]]
    }

    ADDING_DEPTH(
        addingdepthin_ch
    )

    FINAL_REPORT(
        (ADDING_DEPTH.out.repdepth.map{it[1]}).collect()
    )

    BAM_READCOUNT(
        MEDAKA.out.bamfiles.join(REFFIX_FASTA.out.fixedseqref)
    )

    ch_versions = ch_versions.mix(FASTP.out.versions.first())
    ch_versions = ch_versions.mix(NANOPLOT.out.versions.first())
    ch_versions = ch_versions.mix(METAMAPS_MAP.out.versions.first())
    ch_versions = ch_versions.mix(METAMAPS_CLASSIFY.out.versions.first())
    ch_versions = ch_versions.mix(R_METAPLOT.out.versions.first())
    //ch_versions = ch_versions.mix(KRONA_KRONADB.out.versions.first())
    //ch_versions = ch_versions.mix(KRONA_KTIMPORTTAXONOMY.out.versions.first())
    ch_versions = ch_versions.mix(SEQTK_SUBSEQ.out.versions.first())
    ch_versions = ch_versions.mix(MEDAKA.out.versions.first())
    ch_versions = ch_versions.mix(SAMTOOLS_COVERAGE.out.versions.first())
    ch_versions = ch_versions.mix(IVAR_CONSENSUS.out.versions.first())
    ch_versions = ch_versions.mix(HOMOPOLISH_POLISHING.out.versions.first())
    ch_versions = ch_versions.mix(HUMAN_MAPPING.out.versionsmini)
    ch_versions = ch_versions.mix(HUMAN_MAPPING.out.versionssamsort)
    ch_versions = ch_versions.mix(HUMAN_MAPPING.out.versionssamfastq)

    // Wait for RCOVERAGE to complete before running CUSTOM_DUMPSOFTWAREVERSIONS
    CUSTOM_DUMPSOFTWAREVERSIONS(
    ch_versions.unique().collectFile(name: 'collated_versions.yml'),
    ch_rcoverage_done // Add this channel as an input
    )

    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowMetatropics.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowMetatropics.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.json.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(NANOPLOT.out.txt.collect{it[1]}.ifEmpty([]))

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
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
