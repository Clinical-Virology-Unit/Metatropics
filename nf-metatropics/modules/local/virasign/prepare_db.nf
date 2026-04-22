process VIRASIGN_DB {
    tag "Create viral database"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'library://jansendaan94_v2/metatropics/virasign:latest':
        'daanjansen94/virasign:latest' }"

    // Bind DB dir into the container.
    def pipelineRoot = new File("${projectDir}").parentFile.absolutePath
    def db_dir = params.virasign_db_dir ?: "${pipelineRoot}/Databases"
    def dbAbs = file(db_dir).toAbsolutePath().toString()
    if (workflow.containerEngine == 'docker' || workflow.containerEngine == 'podman') {
        containerOptions "-v ${dbAbs}:${dbAbs}"
    }
    if (workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer') {
        containerOptions "--bind ${dbAbs}:${dbAbs}"
    }

    output:
    path "db_ready.txt", emit: ready
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Ensure host path exists.
    file(db_dir).mkdirs()
    def opt = []
    def add = { c, f -> if (c) opt << f }

    // Default DB if unset.
    def rawDbArg = params.virasign_database?.toString()?.trim()
    def effectiveDbArg = rawDbArg ?: 'RVDB'
    add(effectiveDbArg, "-d '${effectiveDbArg}'")
    add(params.virasign_rvdb_version != null, "--rvdb-version ${params.virasign_rvdb_version}")
    add(params.virasign_accessions?.toString()?.trim(), "-a '${params.virasign_accessions}'")
    add(params.virasign_enable_clustering == true, '--enable-clustering')
    add(params.virasign_cluster_identity != null, "--cluster_identity ${params.virasign_cluster_identity}")
    add(params.virasign_max_ambiguous_fraction != null, "--max-ambiguous-fraction ${params.virasign_max_ambiguous_fraction}")

    def tail = task.ext.args?.toString()?.trim()
    def cmd = (['virasign', '--prepare-db', "--db-dir", "${db_dir}"] + opt + (tail ? [tail] : [])).join(' ')

    """
    ${cmd}

    echo "ready" > db_ready.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        virasign: \$(virasign --version 2>&1 | head -n1 || echo 'unknown')
    END_VERSIONS
    """
}

