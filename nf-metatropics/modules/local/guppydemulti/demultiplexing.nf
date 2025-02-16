process GUPPYDEMULTI_DEMULTIPLEXING {
    label 'process_medium'

    container "$projectDir/images/guppy.sif"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        '$projectDir/images/guppy.sif' :
        'daanjansen94/guppy:v6.5.7' }"

    input:
    path reads

    output:
    path("barcoding_summary.txt"), emit: summary
    path("basename.txt"), emit: barcodeList
    path("barcode*"), emit: barcodeReads

    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    if [[ $params.pair == "true" ]]; then
        guppy_barcoder -i $reads -s $params.outdir --require_barcodes_both_ends --num_barcoding_buffers 4 --records_per_fastq 0 -x "auto"

    elif [[ $params.pair == "false" ]]; then
        guppy_barcoder -i $reads -s $params.outdir --num_barcoding_buffers 4 --records_per_fastq 0 -x "auto"
    fi

    cat $params.outdir/barcoding_summary.txt > barcoding_summary.txt
    ls -d $params.outdir/barcode* > list.txt
    for i in `cat list.txt`; do basename \$i ;done > basename.txt
    mv $params.outdir/barcode* .
    rm list.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        guppy_barcoder: \$(echo \$(guppy_barcoder --version) | head -n 1| perl -p -e 's/.+Version //g' | cut -d ' ' -f1)
    END_VERSIONS
    """
}
