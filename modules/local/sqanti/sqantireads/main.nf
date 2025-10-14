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
    tuple val(meta), path("*.html"), emit: html
    path "*_mqc.png"               , emit: multiqc
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
        --annotation bambu_annotation.gtf \\
        --report html \\
        $args

    # extract the second plot to png for multiqc
        grep -o 'data:image/png;base64,[^"]*' sqantiReads_plots.html | sed -n '2p' | sed 's/data:image\\/png;base64,//' | base64 -d > sqanitReads_mqc.png

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

// pca.vals_mqc.tsv

// #id: '_deseq2_pca'
// #section_name: ' DESeq2 PCA plot'
// #description: "PCA plot between samples in the experiment.
// #              These values are calculated using <a href='https://bioconductor.org/packages/release/bioc/html/DESeq2.html'>DESeq2</a>
// #              in the <a href='https://github.com/nf-core/rnaseq/blob/master/bin/deseq2_qc.r'><code>deseq2_qc.r</code></a> script."
// #parent_id: 'sample-relationships'
// #parent_name: 'Sample relationships'
// #parent_description: 'Plots interrogating sample relationships, based on final count matrices'
// #plot_type: 'scatter'
// #anchor: '_deseq2_pca'
// #pconfig:
// #    title: 'DESeq2: Principal component plot'
// #    xlab: PC1
// #    ylab: PC2
// "sample"        "PC1: 50% variance"     "PC2: 15% variance"
// "DE_C_H_15"     -62.9202242220269       -29.4119474775929
// "DE_C_H_16"     -42.7194679133331       25.7169528902837
// "DE_C_H_13"     -39.9764950984932       27.3212767954681
// "DE_C_H_14"     -32.1986742039744       35.8353716214637
// "DE_H_H_19"     -27.7972636041262       -56.8915099345402
// "DE_H_H_20"     65.2732060985042        -11.3737531415994
// "DE_H_H_18"     66.5394688381687        -0.735548805395123
// "DE_H_H_17"     73.7994501052807        9.53915805191234
