process SPACERANGER_COUNT_TO_ANNDATA {
    tag "${meta.id}"
    label "process_low"

    container "quay.io/biocontainers/scanpy:1.7.2"

    input:
    tuple val (meta), path("${meta.id}/*")

    output:
    tuple val(meta), path("st_adata_raw.h5ad"), emit: st_raw
    path("versions.yml"), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    read_st_data.py \\
        --SRCountDir "${meta.id}" \\
        --outAnnData st_adata_raw.h5ad

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scanpy: \$(python -c "import scanpy; print(scanpy.__version__)")
    END_VERSIONS
    """
}