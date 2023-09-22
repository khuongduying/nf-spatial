process UNTAR {
    tag "$archive"
    label 'process_single'

    container "quay.io/nf-core/ubuntu:20.04"

    input:
    tuple val(meta), path(archive)

    output:
    tuple val(meta), path("$prefix"), emit: untar
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args  ?: ''
    def args2 = task.ext.args2 ?: ''
    prefix = task.ext.prefix ?: (meta.id ? "${meta.id}" : archive.baseName.toString().replaceFirst(/\.tar$/, ""))

    """
    mkdir $prefix

    # Check if tar contents have only one top-level directory
    ONE_TOP_DIR=\$(tar -taf ${archive} | grep -o -P "^.*?\\/" | uniq | wc -l)

    # Conditionally untar based on the structure of the archive
    if [[ \$ONE_TOP_DIR -eq 1 ]]; then
        tar -C $prefix --strip-components 1 -xavf $args $archive $args2
    else
        tar -C $prefix -xavf $args $archive $args2
    fi

    # Generate versions file
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        untar: \$(echo \$(tar --version 2>&1) | sed 's/^.*(GNU tar) //; s/ Copyright.*\$//')
    END_VERSIONS
    """
}
