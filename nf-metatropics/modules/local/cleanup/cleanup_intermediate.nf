process CLEANUP_INTERMEDIATE {
    label 'process_single'
    
    input:
    val(dummy_input)

    output:
    path "cleanup_complete.txt", emit: cleanup_done

    script:
    """
    if [ "${workflow.containerEngine}" = "docker" ]; then
        # Remove stopped containers
        docker ps -aq --filter status=exited | xargs -r docker rm -f
        
        # Remove dangling images
        docker images -q --filter dangling=true | xargs -r docker rmi -f
        
        # Remove unused images (keeping those used by running containers)
        docker image prune -af --filter "label!=keep"
        
        echo "Intermediate Docker cleanup completed" > cleanup_complete.txt
    else
        echo "Docker cleanup skipped as Docker is not the container engine" > cleanup_complete.txt
    fi
    """
}
