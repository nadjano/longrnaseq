process MERGE_COUNTS {
    tag "$meta.id"
    label 'process_single'

    // conda "conda-forge::python=3.9 conda-forge::pandas=1.5.3"
    container 'oras://community.wave.seqera.io/library/pandas:2.2.3--e136a7b7218cc69c'

    input:
    tuple val(meta), path(count_files)

    output:
    path("merged_counts.tsv"), emit: merged_counts
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    ${baseDir}/bin/merge_count_files.py \\
        --input_files ${count_files.join(' ')} \\
        --sample_ids ${meta.id.join(' ')} \\
        --output merged_counts.tsv \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """

    stub:
    """
    touch merged_counts.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """
}
