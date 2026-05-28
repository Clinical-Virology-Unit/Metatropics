//
// Check input samplesheet and get read channels
//

include { FIX_NAMES } from '../../../modules/local/fix_names'

workflow FIX {
    take:
    basecalled // file: /path/to/samplesheet.csv

    main:
    FIX_NAMES ( basecalled )
        .map { sample, file -> 
            def meta = [:]
            meta.id = sample
            meta.single_end = true
            tuple(meta, file)
        }
        .set { reads }

    emit:
    reads                                     // channel: [ val(meta), [ reads ] ]
    //versions = FIX_NAMES.out.versions // channel: [ versions.yml ]
}

// Function to get list of [ [id:sample, single_end=true], [file ] ]
//def create_map(sample,fastqin) {
def create_map(sample) {
    def meta = [:]
    meta.id = sample
    meta.single_end = true // Set this to true or false depending on whether the sample is single-end or paired-end
    //def fastq_meta = []
    //if (!file(fastqin).exists()) {
    //    exit 1, "ERROR: Please check file ${fastqin}"
    //}
    //fastq_meta = [ meta , [file(fastqin)] ]
    //meta.path = file(fastqin)
    return meta
    //return[fastq_meta]
}
