//
// Quantification with Oarfish
//

include { OARFISH          } from '../../../modules/local/oarfish/main'
include { CUSTOM_TX2GENE   } from '../../../modules/nf-core/custom/tx2gene'
include { TXIMETA_TXIMPORT } from '../../../modules/nf-core/tximeta/tximport'

include { SUMMARIZEDEXPERIMENT_SUMMARIZEDEXPERIMENT as SE_GENE_UNIFIED       } from '../../../modules/nf-core/summarizedexperiment/summarizedexperiment'
include { SUMMARIZEDEXPERIMENT_SUMMARIZEDEXPERIMENT as SE_TRANSCRIPT_UNIFIED } from '../../../modules/nf-core/summarizedexperiment/summarizedexperiment'

workflow QUANTIFY_PSEUDO_ALIGNMENT {
    take:
    samplesheet               // channel: [ val(meta), /path/to/samplsheet ]
    alignment                 // channel: [ val(meta), [ bam ] ]
    pseudo_aligner            // channel: [ val(meta), 'salmon' | 'kallisto' ]
    gtf                       // channel: [ val(meta), [ gtf ] ]
    gtf_id_attribute          //     val: GTF gene ID attribute
    gtf_extra_attribute       //     val: GTF alternative gene attribute (e.g. gene_name)

    main:
    ch_versions = Channel.empty()

    //
    // Quantify and merge counts across samples
    //
    ch_pseudo_results = OARFISH (
        alignment

    )

    ch_pseudo_results = OARFISH.out.results
    ch_pseudo_multiqc = ch_pseudo_results
    ch_versions = ch_versions.mix(OARFISH.out.versions.first())

    CUSTOM_TX2GENE (
        gtf,
        ch_pseudo_results.collect{ it[1] }.map { [ [:], it ] },
        'salmon',
        gtf_id_attribute,
        gtf_extra_attribute
    )
    ch_versions = ch_versions.mix(CUSTOM_TX2GENE.out.versions)

    TXIMETA_TXIMPORT (
        ch_pseudo_results.collect{ it[1] }.map { [ ['id': 'all_samples'], it ] },
        CUSTOM_TX2GENE.out.tx2gene,
        pseudo_aligner
    )

    ch_versions = ch_versions.mix(TXIMETA_TXIMPORT.out.versions)

    ch_gene_unified = TXIMETA_TXIMPORT.out.counts_gene
                        .join(TXIMETA_TXIMPORT.out.counts_gene_length_scaled)
                        .join(TXIMETA_TXIMPORT.out.counts_gene_scaled)
                        .join(TXIMETA_TXIMPORT.out.lengths_gene)
                        .join(TXIMETA_TXIMPORT.out.tpm_gene)
                        .map{tuple(it[0], it.tail())}

    // SE_GENE_UNIFIED (
    //     ch_gene_unified,
    //     CUSTOM_TX2GENE.out.tx2gene,
    //     samplesheet
    // )
    // ch_versions = ch_versions.mix(SE_GENE_UNIFIED.out.versions)

    ch_transcript_unified = TXIMETA_TXIMPORT.out.counts_transcript
                        .join(TXIMETA_TXIMPORT.out.lengths_transcript)
                        .join(TXIMETA_TXIMPORT.out.tpm_transcript)
                        .map{tuple(it[0], it.tail())}



    // SE_TRANSCRIPT_UNIFIED (
    //     ch_transcript_unified,
    //     CUSTOM_TX2GENE.out.tx2gene,
    //     samplesheet
    // )
    // ch_versions = ch_versions.mix(SE_TRANSCRIPT_UNIFIED.out.versions)

    emit:
    results                       = ch_pseudo_results                              // channel: [ val(meta), results_dir ]
    multiqc                       = ch_pseudo_multiqc                              // channel: [ val(meta), files_for_multiqc ]

    tpm_gene                      = TXIMETA_TXIMPORT.out.tpm_gene                  //    path: *gene_tpm.tsv
    counts_gene                   = TXIMETA_TXIMPORT.out.counts_gene               //    path: *gene_counts.tsv
    lengths_gene                  = TXIMETA_TXIMPORT.out.lengths_gene              //    path: *gene_lengths.tsv
    counts_gene_length_scaled     = TXIMETA_TXIMPORT.out.counts_gene_length_scaled //    path: *gene_counts_length_scaled.tsv
    counts_gene_scaled            = TXIMETA_TXIMPORT.out.counts_gene_scaled        //    path: *gene_counts_scaled.tsv
    tpm_transcript                = TXIMETA_TXIMPORT.out.tpm_transcript            //    path: *gene_tpm.tsv
    counts_transcript             = TXIMETA_TXIMPORT.out.counts_transcript         //    path: *transcript_counts.tsv
    lengths_transcript            = TXIMETA_TXIMPORT.out.lengths_transcript        //    path: *transcript_lengths.tsv

    // merged_gene_rds_unified       = SE_GENE_UNIFIED.out.rds                        //    path: *.rds
    // merged_transcript_rds_unified = SE_TRANSCRIPT_UNIFIED.out.rds                  //    path: *.rds

    versions                      = ch_versions                                    // channel: [ versions.yml ]
}
