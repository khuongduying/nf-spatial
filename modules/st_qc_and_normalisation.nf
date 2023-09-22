process ST_QC_AND_NORMALISATION {
    tag "${meta.id}"
    label 'process_low'

    container "docker.io/erikfas/spatialtranscriptomics"

    input:
    tuple val(meta), path(st_raw, stageAs: "adata_raw.h5ad")

    output:
    tuple val(meta), path("st_adata_norm.h5ad"), emit: st_data_norm
    tuple val(meta), path("st_adata_plain.h5ad"), emit: st_data_plain
    path("versions.yml"), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    st_qc_and_normalisation.py \\
        --raw_adata ${st_raw} \\
        --min_counts ${params.st_preprocess_min_counts} \\
        --min_genes ${params.st_preprocess_min_genes} \\
        --min_cells ${params.st_preprocess_min_cells} \\
        --data_plain_name st_adata_plain.h5ad \\
        --data_norm_name st_adata_norm.h5ad

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scanpy: \$(python -c "import scanpy; print(scanpy.__version__)")
    END_VERSIONS
    """
}