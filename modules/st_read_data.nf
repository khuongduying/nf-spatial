process ST_READ_DATA {
    tag "${meta.id}"
    label "process_low"

    container "docker.io/gcfntnu/scanpy:1.9.2"

    input:
    tuple val (meta), path("${meta.id}/*")

    output:
    tuple val(meta), path("st_adata_raw.h5ad"), emit: st_raw
    path("versions.yml"), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    st_data_read.py \\
        --SR_count_dir "${meta.id}" \\
        --out_AnnData st_adata_raw.h5ad

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scanpy: \$(python -c "import scanpy; print(scanpy.__version__)")
    END_VERSIONS
    """
}