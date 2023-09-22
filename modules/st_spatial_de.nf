// Define the ST_SPATIAL_DE process
process ST_SPATIAL_DE {
    tag "${meta.id}"
    label 'process_medium'

    // Docker container with required software
    container "docker.io/erikfas/spatialtranscriptomics"

    input:
    tuple val(meta), path(st_adata_norm, stageAs: "adata_norm.h5ad")

    output:
    tuple val(meta), path("*.csv"), emit: degs
    path("versions.yml"), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    st_spatial_de.py \\
        --file_name_st ${st_adata_norm} \\
        --save_spatial_de_file_name st_spatial_de.csv

    // Generate a YAML file containing version information of used software
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        leidenalg: \$(python -c "import leidenalg; print(leidenalg.version)")
        scanpy: \$(python -c "import scanpy; print(scanpy.__version__)")
        SpatialDE: \$(python -c "from importlib.metadata import version; print(version('SpatialDE'))")
    END_VERSIONS
    """
}