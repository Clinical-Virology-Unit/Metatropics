//
// Check input samplesheet and get read channels
//

include { MINIMAP2_ALIGN              } from '../../../modules/nf-core/minimap2/align/main'
include { SAMTOOLS_SORT               } from '../../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_hFASTQ            } from '../../../modules/nf-core/samtools/fastq/main'

workflow HUMAN_MAPPING {
    take:
    readsONT
    host_fasta

    main:
    MINIMAP2_ALIGN(
        readsONT,
        host_fasta,
        true,
        false,
        false
    )
    SAMTOOLS_SORT(
        MINIMAP2_ALIGN.out.bam
    )

    SAMTOOLS_hFASTQ(
        SAMTOOLS_SORT.out.bam,
        false
    )

    emit:
    humanout = SAMTOOLS_hFASTQ.out.other  // Changed from nohumanreads to humanout
    versionsmini = MINIMAP2_ALIGN.out.versions
    versionssamsort = SAMTOOLS_SORT.out.versions
    versionssamfastq = SAMTOOLS_hFASTQ.out.versions
}
