// host_mapping.nf

include { MINIMAP2_ALIGN } from '../../../modules/nf-core/minimap2/align/main'
include { SAMTOOLS_SORT  } from '../../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_hoFASTQ } from '../../../modules/nf-core/samtools/fastq/main2'

workflow HOST_MAPPING {
    take:
    reads
    host_fasta

    main:
    MINIMAP2_ALIGN(
        reads,
        host_fasta,
        true,
        false,
        false
    )
    SAMTOOLS_SORT(
        MINIMAP2_ALIGN.out.bam
    )

    SAMTOOLS_hoFASTQ(
        SAMTOOLS_SORT.out.bam,
        false
    )

    emit:
    hostout = SAMTOOLS_hoFASTQ.out.other  // Changed from nohostReads to hostout
    versionsMini = MINIMAP2_ALIGN.out.versions
    versionsSamSort = SAMTOOLS_SORT.out.versions
    versionsSamFastq = SAMTOOLS_hoFASTQ.out.versions
}
