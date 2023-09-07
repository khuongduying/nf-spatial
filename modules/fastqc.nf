process FASTQC {
    tag "${meta.id}"
    label 'process_medium'

    container "quay.io/biocontainers/fastqc:0.12.1"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.zip"), emit: zip
    path("versions.yml"), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    // Generate pairs of old and new filenames for renaming before running fastqc
    def filename_pairs = reads instanceof Path || reads.size() == 1 
        ? [[ reads, "${prefix}.${reads.extension}" ]] 
        : reads.withIndex().collect { entry, index -> [ entry, "${prefix}_${index + 1}.${entry.extension}" ] }
    
    def rename_commands = filename_pairs.collect { oldName, newName -> "ln -s \"${oldName}\" \"${newName}\"" }.join('; ')
    def renamed_files = filename_pairs.collect { _, newName -> newName }.join(' ')
    """
    set -e
    $rename_commands
    fastqc $args --threads $task.cpus $renamed_files

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: \$( fastqc --version | sed -e "s/FastQC v//g" )
    END_VERSIONS
    """
}