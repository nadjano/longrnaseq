
process HAPLOTAG {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'community.wave.seqera.io/library/pip_argparse_pathlib_pysam:bdbdbf01755d4994' :
        'community.wave.seqera.io/library/pip_argparse_pathlib_pysam:bdbdbf01755d4994' }"


    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.haplotagged.bam"), emit: bam
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    python \\
        ${projectDir}/bin/primary_alignment_haplotag.py \\
        -i ${bam} \\
        -o ${prefix}.haplotagged.bam \\
        -t HP \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pysam: \$(python -c "import pysam; print(pysam.__version__)")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.haplotagged.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pysam: \$(python -c "import pysam; print(pysam.__version__)")
    END_VERSIONS
    """
}
