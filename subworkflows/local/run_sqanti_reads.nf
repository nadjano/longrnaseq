// TODO nf-core: If in doubt look at other nf-core/subworkflows to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/subworkflows
//               You can also ask for help via your pull request or on the #subworkflows channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A subworkflow SHOULD import at least two modules


include { SPLICEDBAM2GFF                               } from '../../modules/local/splicedbam2gff'
include { SQANTIREADS } from '../../../modules/modules/nf-core/sqantireads'

workflow RUN_SQANTI_READS {

    take:
    // TODO nf-core: edit input (take) channels
    ch_bam // channel: [ val(meta), [ bam ] ]
    ch_fasta // channel: [ val(meta), [ fasta ] ]
    ch_gtf // channel: [ val(meta), [ gtf ] ]
    main:

    ch_versions = Channel.empty()


    //
    // MODULES:
    //
    SPLICEDBAM2GFF  ( ch_bam )
    ch_versions = ch_versions.mix(SPLICEDBAM2GFF.out.versions.first())

    //
    // MODULES: RUN SQANTI READS
    //
    SQANTIREADS (
                SPLICEDBAM2GFF.out.gff,
                ch_fasta,
                ch_gtf
                )


    emit:
    // TODO nf-core: edit emitted channels
    gff      = SPLICEDBAM2GFF.out.gff          // channel: [ val(meta), [ gff ] ]


    versions = ch_versions                     // channel: [ versions.yml ]
}

