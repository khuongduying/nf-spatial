//
// Run Space Ranger count
//

include { UNTAR as SPACERANGER_UNTAR_REFERENCE } from "../modules/untar.nf"
include { SPACERANGER_COUNT } from '../modules/spaceranger_count.nf'

workflow SPACERANGER {
    
    // Input channel for spatial transcriptomics data
    take:
    ch_st_data // channel: [ val(meta), [ raw st data ] ]

    main:
    // Versions tracking
    def ch_versions = Channel.empty()

    // Handling reference files
    def ch_reference = Channel.empty()
    if (params.spaceranger_reference ==~ /.*\.tar\.gz$/) {
        ref_file = file(params.spaceranger_reference)
        SPACERANGER_UNTAR_REFERENCE([[id: "reference"], ref_file])
        ch_reference = SPACERANGER_UNTAR_REFERENCE.out.untar.map({ meta, ref -> ref })
        ch_versions.mix(SPACERANGER_UNTAR_REFERENCE.out.versions)
    } else {
        ch_reference = file(params.spaceranger_reference, type: "dir", checkIfExists: true)
    }

    // Handling optional probeset
    ch_probeset = params.spaceranger_probeset ? 
                  file(params.spaceranger_probeset, checkIfExists: true) : 
                  Channel.empty()
    
    println(ch_st_data)
    println(ch_reference)
    println(ch_probeset)

    // Running Space Ranger count
    SPACERANGER_COUNT(ch_st_data, ch_reference, ch_probeset)
    ch_versions.mix(SPACERANGER_COUNT.out.versions.first())

    // Outputs
    emit:
    sr_dir   = SPACERANGER_COUNT.out.outs
    versions = ch_versions  // channel: [ versions.yml ]
}