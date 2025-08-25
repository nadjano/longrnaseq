// TODO nf-core: If in doubt look at other nf-core/subworkflows to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/subworkflows
//               You can also ask for help via your pull request or on the #subworkflows channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A subworkflow SHOULD import at least two modules


include { SPLICEDBAM2GFF                               } from '../../modules/local/splicedbam2gff'
include { SQANTIREADS                                  } from '../../modules/local/sqanti/sqantireads'
include { SQANTIQC                                     } from '../../modules/local/sqanti/sqantiqc'
include { SAMTOOLS_VIEW as SAMTOOLS_FILTER            } from '../../modules/nf-core/samtools/view/main'
workflow RUN_SQANTI_READS {

    take:
    // TODO nf-core: edit input (take) channels
    ch_bam // channel: [ val(meta), [ bam ], [ bai] ]
    ch_fasta // channel: [ val(meta), [ fasta ] ]
    ch_gtf // channel: [ val(meta), [ gtf ] ]
    main:

    ch_versions = Channel.empty()

    //
    // MODULES: RUN SAMTOOLS VIEW to filter BAM files and downsample
    //
    SAMTOOLS_FILTER (
                ch_bam,
                ch_fasta,
                [], // qname
                'bai' // index_format
                )


    // //
    // // MODULES: RUN SPLICEDBAM2GFF
    // //
    // SPLICEDBAM2GFF  ( SAMTOOLS_FILTER.out.bam
    //                 )
    // ch_versions = ch_versions.mix(SPLICEDBAM2GFF.out.versions.first())


    // //
    // // MODULES: RUN SQANTI QC
    // //
    // // For testing only select two samples
    // // Instead of reassigning, create a conditional channel
    // gff_channel = params.sqanti_test ?
    //     SPLICEDBAM2GFF.out.gff.take(2) :
    //     SPLICEDBAM2GFF.out.gff


    // SQANTIQC (
    //             gff_channel,
    //             ch_fasta,
    //             ch_gtf
    //             )

    // // Collect the gff files
    // combined_sqanti_ch = SQANTIQC.out.sqanti_qc
    // .toSortedList { a, b ->
    //     // Sort by sample ID
    //     a[0].id <=> b[0].id
    // }
    // .map { tuples ->
    //     def metas = tuples.collect { it[0] }  // Extract all meta maps
    //     def files = tuples.collect { it[1] }  // Extract all quant files
    //     [metas, files]
    // }


    // // MODULES: RUN SQANTI READS

    // SQANTIREADS (
    //             combined_sqanti_ch,
    //             ch_gtf
    //             )


    emit:
    // TODO nf-core: edit emitted channels
    // gff      = SPLICEDBAM2GFF.out.gff          // channel: [ val(meta), [ gff ] ]
    bam      = SAMTOOLS_FILTER.out.bam          // channel: [ val(meta), [ bam ], [ bai ] ]

    versions = ch_versions                     // channel: [ versions.yml ]
}

