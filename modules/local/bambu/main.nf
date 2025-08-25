process BAMBU {
    label 'process_high'
    label 'high_memory'

    conda "conda-forge::r-base=4.0.3 bioconda::bioconductor-bambu=3.0.8 bioconda::bioconductor-bsgenome=1.66.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'community.wave.seqera.io/library/bioconductor-bambu_r-xgboost:c4dc8df46f69eaa4' :
        'community.wave.seqera.io/library/bioconductor-bambu_r-xgboost:c4dc8df46f69eaa4' }"



    input:
    path(fasta)
    tuple val(meta), path(gtf)
    path bams

    output:
    path "extended_annotations.gtf", emit: extended_gtf
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    run_bambu.r \\
        --tag=. \\
        --ncore=$task.cpus \\
        --annotation=$gtf \\
        --fasta=$fasta $bams

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bioconductor-bambu: \$(Rscript -e "library(bambu); cat(as.character(packageVersion('bambu')))")
        bioconductor-bsgenome: \$(Rscript -e "library(BSgenome); cat(as.character(packageVersion('BSgenome')))")
        conda-forge-r-xgboost: \$(Rscript -e "library(xgboost); cat(as.character(packageVersion('xgboost')))")
    END_VERSIONS
    """
}
