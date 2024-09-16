process CLEANUP {
    label 'process_single'
    
    input:
    path versions_file
    path final_report
    path read_counts_csv
    
    output:
    path "cleanup_complete.txt", emit: cleanup_done
    
    script:
    """
    if [ "${workflow.containerEngine}" = "docker" ]; then
        docker ps -aq | xargs -r docker rm -f
        docker images -q | xargs -r docker rmi -f
        echo "Final Docker cleanup completed" > cleanup_complete.txt
    else
        echo "Docker cleanup skipped as Docker is not the container engine" > cleanup_complete.txt
    fi
    """
}
