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
                'csi' // index_format, for large genomes csi is required
                )

  
        //
        // MODULES: RUN SPLICEDBAM2GFF
        //
        SPLICEDBAM2GFF  ( SAMTOOLS_FILTER.out.bam
                        )
        ch_versions = ch_versions.mix(SPLICEDBAM2GFF.out.versions.first())


        //
        // MODULES: RUN SQANTI QC
        //

        gff_channel = SPLICEDBAM2GFF.out.gff


        SQANTIQC (
                    gff_channel,
                    ch_fasta,
                    ch_gtf
                    )
        ch_versions = ch_versions.mix(SQANTIQC.out.versions.first())
        // Collect the gff files
        combined_sqanti_ch = SQANTIQC.out.sqanti_qc
        .toSortedList { a, b ->
            // Sort by sample ID
            a[0].id <=> b[0].id
        }
        .map { tuples ->
            def metas = tuples.collect { it[0] }  // Extract all meta maps
            def files = tuples.collect { it[1] }  // Extract all quant files
            [metas, files]
        }


        // MODULES: RUN SQANTI READS

        SQANTIREADS (
                    combined_sqanti_ch,
                    ch_gtf
                    )
        ch_versions = ch_versions.mix(SQANTIREADS.out.versions.first())

        emit:
        // TODO nf-core: edit emitted channels
        gff       = SPLICEDBAM2GFF.out.gff          // channel: [ val(meta), [ gff ] ]
        sqanti_qc = SQANTIQC.out.sqanti_qc         // channel: [ val(meta), [ qc ] ]
        bam       = SAMTOOLS_FILTER.out.bam          // channel: [ val(meta), [ bam ], [ bai ] ]
        multiqc   =  SQANTIREADS.out.multiqc   // channel: [ _mqc.png ]
        versions  = ch_versions                     // channel: [ versions.yml ]

    



}

