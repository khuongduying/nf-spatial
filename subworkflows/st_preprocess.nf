//
// Pre-processing of Spatial Transcriptomics (ST) data
//

// Import the ST_QC_AND_NORMALISATION module
include { ST_QC_AND_NORMALISATION } from '../../modules/local/st_qc_and_normalisation'

// Define the ST_PREPROCESS workflow
workflow ST_PREPROCESS {

    // Input channel
    take:
    st_raw  // Raw spatial transcriptomics data

    main:

    // Initialize channel for software versions used in the pipeline
    ch_versions = Channel.empty()

    //
    // Spatial data pre-processing
    //
    // Run the QC and normalization steps on the raw spatial transcriptomics data.
    ST_QC_AND_NORMALISATION (
        st_raw
    )

    // Collect version information from the QC and Normalisation step
    ch_versions = ch_versions.mix(ST_QC_AND_NORMALISATION.out.versions)

    // Output channels
    emit:
    st_data_norm  = ST_QC_AND_NORMALISATION.out.st_data_norm  // Normalized data
    st_data_plain = ST_QC_AND_NORMALISATION.out.st_data_plain // Data without normalization
    st_html       = ST_QC_AND_NORMALISATION.out.html          // QC and normalization HTML report
    versions      = ch_versions                               // Software versions used

} // End of ST_PREPROCESS workflow