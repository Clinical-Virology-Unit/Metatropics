process FINAL_REPORT {
    //tag "${meta.id}"
    
    input:
    path(report)
    
    output:
    path("all.final_report.tsv"), emit: finalReport
    
    script:
    def horizontal_cov_threshold = params.horizontal_coverage != null ? params.horizontal_coverage : 1
    """
    echo -e "Sample\tAccession number\tVirus\tReads\tVertical coverage\tHorizontal coverage\tRead identity\tRead length\tBase quality" > all.final_report.tsv
    cat *.sdepth.tsv | grep -v VirusName | awk -F'\\t' '{
        gsub(/_T1/, "", \$1)
        gsub(/"/, "", \$1)
        gsub(/"/, "", \$2)
        gsub(/"/, "", \$4)
        if (\$10 > ${horizontal_cov_threshold} || ${horizontal_cov_threshold} == 0) {  # Filter on horizontal coverage > threshold (or include all if threshold is 0)
            print \$1"\\t"\$2"\\t"\$4"\\t"\$5"\\t"\$9"\\t"\$10"\\t"\$12"\\t"\$13"\\t"\$14
        }
    }' >> all.final_report.tsv
    """
}
