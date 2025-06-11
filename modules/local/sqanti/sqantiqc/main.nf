
process SQANTIQC {
    tag "$meta.id"

    cache 'lenient'

    // container 'docker://anaconesalab/sqanti3'
    //container "/scratch/nadjafn/nf-core-plantlongrnaseq/work/singularity/anaconesalab-sqanti3-latest.img"
    // container  'docker://pippel/squanti3:v5.2.1'
    // container "docker://anaconesalab/sqanti3:5.3.6-conda-fix"

    input:
    tuple val(meta), path(gff)
    tuple val(meta2), path(genome)
    path(annotation)

    output:
    tuple val(meta), path("${meta}/*.html"), emit: html
    tuple val(meta), path("${meta}", type: 'dir'), emit: sqanti_qc
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    python ${params.sqanti_dir}/sqanti3_qc.py \\
        --isoforms $gff \\
        --refGTF $annotation \\
        --refFasta $genome \\
        --skipORF --min_ref_len 0 --chunks 10 \\
        -t $task.cpus -d ${prefix} -o ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sqantireads: \$(python ${params.sqanti_dir}/sqanti3_qc.py --version |& sed '1!d ; s/sqanti-reads //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // TODO nf-core: A stub section should mimic the execution of the original module as best as possible
    //               Have a look at the following examples:
    //               Simple example: https://github.com/nf-core/modules/blob/818474a292b4860ae8ff88e149fbcda68814114d/modules/nf-core/bcftools/annotate/main.nf#L47-L63
    //               Complex example: https://github.com/nf-core/modules/blob/818474a292b4860ae8ff88e149fbcda68814114d/modules/nf-core/bedtools/split/main.nf#L38-L54
    """
    touch ${prefix}
    touch ${prefix}/.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sqantireads: \$(python ${params.sqanti_dir}/sqanti3_qc.py --version  |& sed '1!d ; s/sqanti-reads //')
    END_VERSIONS
    """
}
