def checkPathParamList = [ params.input, params.spaceranger_reference, params.spaceranger_probeset ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

include { SPACERANGER_COUNT_TO_ANNDATA } from '../modules/st_read_data'
include { INPUT_CHECK } from '../subworkflows/input_check'
include { SPACERANGER } from '../subworkflows/spaceranger'
include { ST_PREPROCESS } from '../subworkflows/st_preprocess'
include { ST_POSTPROCESS } from '../subworkflows/st_postprocess'
include { FASTQC } from "../modules/fastqc"

// Spatial transcriptomics workflow

workflow ST {
    ch_versions = Channel.empty()

    // SUBWORKFLOW: Read and validate samplesheet
    INPUT_CHECK (
        ch_input
    )

    // MODULE: FastQC
    FASTQC(
        INPUT_CHECK.out.ch_spaceranger_input.map{ it -> [it[0] /* meta */, it[1] /* reads */]}
    )
    ch_versions = ch_versions.mix(FASTQC.out.versions)

    // SUBWORKFLOW: Space Ranger raw data processing
    SPACERANGER (
        INPUT_CHECK.out.ch_spaceranger_input
    )
    ch_versions = ch_versions.mix(SPACERANGER.out.versions)
    ch_downstream_input = INPUT_CHECK.out.ch_downstream_input.concat(SPACERANGER.out.sr_dir).map {
        meta, outs -> [meta, outs.findAll { it -> 
            Utils.DOWNSTREAM_REQUIRED_SPACERANGER_FILES.contains(it.name) 
        }]
    }

    // MODULE: Read ST data and save as `anndata`
    ST_READ_DATA (
        ch_downstream_input
    )
    ch_versions = ch_versions.mix(ST_READ_DATA.out.versions)

    // SUBWORKFLOW: Pre-processing of ST  data
    ST_PREPROCESS (
        ST_READ_DATA.out.st_raw
    )
    ch_versions = ch_versions.mix(ST_PREPROCESS.out.versions)

    // SUBWORKFLOW: Post-processing and reporting
    ST_POSTPROCESS (
        ST_PREPROCESS.out.st_data_norm
    )
    ch_versions = ch_versions.mix(ST_POSTPROCESS.out.versions)

    // MODULE: Pipeline reporting
}