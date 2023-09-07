process SPACERANGER_COUNT {
    tag "${meta.id}"
    label "process_high"

    container "quay.io/cumulus/spaceranger:2.0.1"

    input:
    tuple val(meta), path(reads), path(image), path(cytaimage), path(darkimage), path(colorizedimage), path(alignment), path(slidefile)
    path(reference)
    path(probeset)

    output:
    tuple val(meta), path("**/outs/**"), emit: outs
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    // Add flags for optional inputs on demand.
    def probeset = probeset ? "--probe-set=\"${probeset}\"" : ""
    def alignment = alignment ? "--loupe-alignment=\"${alignment}\"" : ""
    def slidefile = slidefile ? "--slidefile=\"${slidefile}\"" : ""
    def image = image ? "--image=\"${image}\"" : ""
    def cytaimage = cytaimage ? "--cytaimage=\"${cytaimage}\"" : ""
    def darkimage = darkimage ? "--darkimage=\"${darkimage}\"" : ""
    def colorizedimage = colorizedimage ? "--colorizedimage=\"${colorizedimage}\"" : ""
    """
    spaceranger count \\
        --id="${prefix}" \\
        --sample="${meta.id}" \\
        --fastqs=. \\
        --slide="${meta.slide}" \\
        --area="${meta.area}" \\
        --transcriptome="${reference}" \\
        --localcores=${task.cpus} \\
        --localmem=${task.memory.toGiga()} \\
        $image $cytaimage $darkimage $colorizedimage \\
        $probeset \\
        $alignment \\
        $slidefile \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        spaceranger: \$(spaceranger -V | sed -e "s/spaceranger spaceranger-//g")
    END_VERSIONS
    """
}