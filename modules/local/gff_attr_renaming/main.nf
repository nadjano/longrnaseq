process GFF_ATTR_RENAMING {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pyranges:0.1.2--pyhdfd78af_1' :
        'biocontainers/pyranges:0.1.4--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(gtf)


    output:
    tuple val(meta), path("*.fixed.gtf"), emit: fixed_gtf
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args  = task.ext.args  ?: '' // args is used for the main arguments of the tool
    def args2 = task.ext.args2 ?: '' // args2 is used to specify a program when no program file has been given
    prefix    = task.ext.prefix ?: "${meta.id}"

    """
    python ${projectDir}/bin/fix_gtf.py $gtf ${prefix}.fixed.gtf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gff_attr_renaming: \$(python --version | sed 's/Python //')
        pyranges: \$(python -c "import pyranges; print(pyranges.__version__)")
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"



    """
    touch ${prefix}.fixed.gtf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gff_attr_renaming: \$(python --version | sed 's/Python //')
        pyranges: \$(python -c "import pyranges; print(pyranges.__version__)")
    END_VERSIONS
    """
}
