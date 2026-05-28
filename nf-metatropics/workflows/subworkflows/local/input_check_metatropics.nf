//
// Check input samplesheet and get read channels
//

include { SAMPLESHEET_CHECK_METATROPICS } from '../../../modules/local/samplesheet_check_metatropics'

workflow INPUT_CHECK_METATROPICS {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    SAMPLESHEET_CHECK_METATROPICS ( samplesheet )
        .csv
        .splitCsv ( header:true, sep:',' )
        .map { create_fastq_channel(it) }
        .set { reads }

    emit:
    reads                                     // channel: [ val(meta), [ reads ] ]
    versions = SAMPLESHEET_CHECK_METATROPICS.out.versions // channel: [ versions.yml ]
}

// Function to get list of [ sample_id, barcode_path_or_label ]
def create_fastq_channel(LinkedHashMap row) {
    def meta = [:]
    meta.id         = row.sample
    meta.single_end = true  // Nanopore: always single-end; paired-end is not supported.

    return [ meta.id, row.barcode ]
}
