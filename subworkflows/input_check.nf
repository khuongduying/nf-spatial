// Function to validate and create an input channel for downstream processing
def create_channel_downstream(LinkedHashMap meta) {
    meta["id"] = meta.remove("sample")  // Renaming 'sample' to 'id'
    spaceranger_dir = file("${meta.remove('spaceranger_dir')}/**")  // Directory with required spaceranger files

    DOWNSTREAM_REQUIRED_SPACERANGER_FILES = [
        "raw_feature_bc_matrix.h5",
        "scalefactors_json.json",
        "tissue_hires_image.png",
        "tissue_lowres_image.png"
    ]

    // Check for the tissue_positions file. 
    // Either "tissue_positions.csv" or "tissue_positions_list.csv" is acceptable.
    tissue_positions_present = spaceranger_dir*.name.contains("tissue_positions.csv") || spaceranger_dir*.name.contains("tissue_positions_list.csv")

    if (!tissue_positions_present) {
        error "The specified spaceranger output directory doesn't contain either 'tissue_positions.csv' or 'tissue_positions_list.csv' for sample `${meta.id}`"
    }

    // Loop over the other expected spaceranger files and check their presence
    for (f in DOWNSTREAM_REQUIRED_SPACERANGER_FILES) {
        if (!spaceranger_dir*.name.contains(f)) {
            error "The specified spaceranger output directory doesn't contain the required file `${f}` for sample `${meta.id}`"
        }
    }

    // Return metadata and directory path
    return [meta, spaceranger_dir]
}

// Function to validate and create an input channel for spaceranger processing
def create_channel_spaceranger(LinkedHashMap meta) {
    meta["id"] = meta.remove("sample")  // Renaming 'sample' to 'id'
    
    // Utility function to fetch a file path from meta and check its existence
    def get_file_from_meta = { key ->
        v = meta.remove(key);
        return v ? file(v) : []
    }

    // Extract and check fastq files
    fastq_dir = meta.remove("fastq_dir")
    fastq_files = files(fastq_dir + "/*${meta['id']}*.fastq.gz")
    if (!fastq_files.size()) {
        error "No `fastq_dir` specified or no samples found in folder."
    } else {
        log.info "${fastq_files.size()} FASTQ files found for sample ${meta['id']}."
    }

    // Extract and validate optional files
    manual_alignment = get_file_from_meta("manual_alignment")
    slidefile = get_file_from_meta("slidefile")
    image = get_file_from_meta("image")
    cytaimage = get_file_from_meta("cytaimage")
    colorizedimage = get_file_from_meta("colorizedimage")
    darkimage = get_file_from_meta("darkimage")

    check_optional_files = ["manual_alignment", "slidefile", "image", "cytaimage", "colorizedimage", "darkimage"]
    for (k in check_optional_files) {
        if (this.binding[k] && !this.binding[k].exists()) {
            error "File for `${k}` is specified, but does not exist: ${this.binding[k]}."
        }
    }
    if (!(image || cytaimage || colorizedimage || darkimage)) {
        error "Need to specify at least one of 'image', 'cytaimage', 'colorizedimage', or 'darkimage' in the samplesheet"
    }

    // Return all validated inputs
    return [meta, fastq_files, image, cytaimage, darkimage, colorizedimage, manual_alignment, slidefile]
}

workflow INPUT_CHECK {
    // Inputs for the workflow
    take:
    samplesheet

    main:
    // Core logic
    // Reads the samplesheet and splits the CSV by rows. It then categorizes each row
    // into one of two channels based on the presence of 'spaceranger_dir'
    ch_st = Channel.from(samplesheet).splitCsv(
        header: true,
        sep: ','
    ).branch {
        spaceranger: !it.containsKey("spaceranger_dir")  // Entries without 'spaceranger_dir'
        downstream: it.containsKey("spaceranger_dir")  // Entries with 'spaceranger_dir'
    }

    // Maps each entry from the respective channels into two methods
    ch_spaceranger_input = ch_st.spaceranger.map { create_channel_spaceranger(it) }
    ch_downstream_input = ch_st.downstream.map { create_channel_downstream(it) }

    // Emit the outputs from this workflow for further use
    emit:
    ch_spaceranger_input
    ch_downstream_input
}
