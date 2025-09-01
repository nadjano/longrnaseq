process OARFISH {
    tag "$meta.id"
    label 'process_single'

    cache 'lenient'

    // TODO nf-core: List required Conda package(s).
    //               Software MUST be pinned to channel (i.e. "bioconda"), version (i.e. "1.10").
    //               For Conda, the build (i.e. "h9402c20_2") must be EXCLUDED to support installation on different operating systems.
    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'community.wave.seqera.io/library/oarfish:0.9.0--4958b2b67f1dd777':
        'biocontainers/oarfish:0.8.1--h5ca1c30_0' }"


    input:
    // TODO nf-core: Where applicable all sample-specific information e.g. "id", "single_end", "read_group"
    //               MUST be provided as an input via a Groovy Map called "meta".
    //               This information may not be required in some instances e.g. indexing reference genome files:
    //               https://github.com/nf-core/modules/blob/master/modules/nf-core/bwa/index/main.nf
    // TODO nf-core: Where applicable please provide/convert compressed files as input/output
    //               e.g. "*.fastq.gz" and NOT "*.fastq", "*.bam" and NOT "*.sam" etc.
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${prefix}")        , emit: results
    tuple val(meta), path("${meta.id}.ambig_info.tsv")        , emit: ambig_info
    tuple val(meta), path("*.meta_info.json")  , emit: json_info
    tuple val(meta), path(".command.log")            , emit: log
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir $prefix && oarfish $args \\
            -j $task.cpus \\
            --alignments $bam \\
            -o $prefix \\
            --filter-group no-filters \\
            --write-assignment-probs \\
            --score-threshold 1.0

    # Add the transcript ids to the ambig_info.tsv file
    paste <(awk '{print \$1}' "${meta.id}.quant") \\
      <(awk '{print \$1 "\\t" \$2}' "${meta.id}.ambig_info.tsv") \\
      > "${meta.id}.ambig_info.tmp" && \\
    mv "${meta.id}.ambig_info.tmp" "${meta.id}.ambig_info.tsv"

    # Rename .qant file to qant.sf to match salmon output
    mv ${prefix}.quant ${prefix}/quant.sf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        oarfish: \$(oarfish --version |& sed '1!d ; s/oarfish //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    // TODO nf-core: A stub section should mimic the execution of the original module as best as possible
    //               Have a look at the following examples:
    //               Simple example: https://github.com/nf-core/modules/blob/818474a292b4860ae8ff88e149fbcda68814114d/modules/nf-core/bcftools/annotate/main.nf#L47-L63
    //               Complex example: https://github.com/nf-core/modules/blob/818474a292b4860ae8ff88e149fbcda68814114d/modules/nf-core/bedtools/split/main.nf#L38-L54
    """
    mkdir -p $prefix
    touch ${prefix}.log
    touch ${prefix}.meta_info.json
    touch ${prefix}.ambig_info.json


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        oarfish: \$(oarfish --version |& sed '1!d ; s/oarfish //')
    END_VERSIONS
    """
}
