process SQANTIREADS {

    tag "$meta.id"
    label 'process_medium'


    //container 'docker://anaconesalab/sqanti3'

    input:
    // TODO nf-core: Where applicable all sample-specific information e.g. "id", "single_end", "read_group"
    //               MUST be provided as an input via a Groovy Map called "meta".
    //               This information may not be required in some instances e.g. indexing reference genome files:
    //               https://github.com/nf-core/modules/blob/master/modules/nf-core/bwa/index/main.nf
    // TODO nf-core: Where applicable please provide/convert compressed files as input/output
    //               e.g. "*.fastq.gz" and NOT "*.fastq", "*.bam" and NOT "*.sam" etc.
    tuple val(meta), path(sqanti_dir)
    path(annotation)

    output:
    // TODO nf-core: Named file extensions MUST be emitted for ALL output channels
    tuple val(meta), path("*.pdf"), emit: pdf
    // TODO nf-core: List additional required output channels/values here
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    # make design file
    echo -e "sampleID,file_acc" > design.csv

    # Process each gff file
    for gff_file in ${sqanti_dir}; do
        fbname=\$(basename \$gff_file .gff)
        echo -e "\$fbname,\$fbname" >> design.csv
    done

    cat design.csv

    python ${params.sqanti_dir}/sqanti3_reads.py \\
        --design design.csv \\
        --annotation $annotation \\
        $args \\
        --force_id_ignore

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sqantireads: \$(python ${params.sqanti_dir}/sqanti3_reads.py --version |& sed '1!d ; s/sqanti-reads //')
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
    touch ${prefix}.html


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sqantireads: \$(python ${params.sqanti_dir}/sqanti3_reads.py --version  |& sed '1!d ; s/sqanti-reads //')
    END_VERSIONS
    """
}
