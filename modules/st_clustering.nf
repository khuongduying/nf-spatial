process ST_CLUSTERING {
    tag "${meta.id}"
    label 'process_low'

    container "docker.io/erikfas/spatialtranscriptomics"

    input:
    tuple val(meta), path(st_adata_norm, stageAs: "adata_norm.h5ad")

    output:
    tuple val(meta), path("st_adata_processed.h5ad"), emit: st_adata_processed
    path("versions.yml"), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    st_clustering.py \\
        --file_name_st ${st_adata_norm} \\
        --resolution ${params.st_cluster_resolution} \\
        --save_file_st st_adata_processed.h5ad

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scanpy: \$(python -c "import scanpy; print(scanpy.__version__)")
    END_VERSIONS
    """
}