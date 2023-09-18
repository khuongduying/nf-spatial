process ST_QC_AND_NORMALISATION {
    tag "${meta.id}"
    label 'process_low'

    container "docker.io/erikfas/spatialtranscriptomics"

    input:
    tuple val(meta), path(st_raw, stageAs: "adata_raw.h5ad")

    output:
    tuple val(meta), path("st_adata_norm.h5ad"), emit: st_data_norm
    tuple val(meta), path("st_adata_plain.h5ad"), emit: st_data_plain
    tuple val(meta), path("st_qc_and_normalisation.html"), emit: html
    path("versions.yml"), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    st_qc_and_normalisation.py \\
        --raw-adata ${st_raw} \\
        --min-counts ${params.st_preprocess_min_counts} \\
        --min-genes ${params.st_preprocess_min_genes} \\
        --min-cells ${params.st_preprocess_min_cells} \\
        --data-plain-name st_adata_plain.h5ad \\
        --data-norm-name st_adata_norm.h5ad

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scanpy: \$(python -c "import scanpy; print(scanpy.__version__)")
    END_VERSIONS
    """
}