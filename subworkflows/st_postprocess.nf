//
// Post-processing of Spatial Transcriptomics (ST) data
// Includes spatial differential expression and clustering
//

// Import required modules
include { ST_SPATIAL_DE } from '../modules/st_spatial_de'
include { ST_CLUSTERING } from '../modules/st_clustering'

// Define the ST_POSTPROCESS workflow
workflow ST_POSTPROCESS {

    // Input channel
    take:
    st_adata_norm  // Normalized spatial transcriptomics data

    main:

    // Initialize channel for collecting software versions used in the pipeline
    ch_versions = Channel.empty()

    //
    // Clustering
    //
    // Perform clustering on the normalized spatial transcriptomics data.
    // The report and normalized data are used as input parameters.
    //
    ST_CLUSTERING (
        st_adata_norm
    )
    // Collect version information from the clustering step
    ch_versions = ch_versions.mix(ST_CLUSTERING.out.versions)

    //
    // Spatial Differential Expression
    //
    // Perform spatial differential expression analysis.
    // It takes as input the report and processed data from the clustering step.
    //
    ST_SPATIAL_DE (
        ST_CLUSTERING.out.st_adata_processed
    )
    // Collect version information from the spatial differential expression step
    ch_versions = ch_versions.mix(ST_SPATIAL_DE.out.versions)

    // Output channels
    emit:
    spatial_degs = ST_SPATIAL_DE.out.degs  // Differential expression genes (DEGs)
    versions     = ch_versions             // Software versions used in the pipeline

} // End of ST_POSTPROCESS workflow