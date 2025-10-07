//
// Run SAMtools stats, flagstat and idxstats
//
// Run SAMtools stats, flagstat and idxstats
//

include { HAPLOTAG                                        } from '../../modules/local/haplotag/main'
include { SAMTOOLS_VIEW as SAMTOOLS_DOWNSAMPLE            } from '../../modules/nf-core/samtools/view/main'

workflow ALIGNMENT_VISUALISATION {
    take:
    ch_bam     // channel: [ val(meta), path(bam), path(bai) ]
    ch_fasta   // channel: [ val(meta), path(fasta) ]

    main:
    ch_versions = Channel.empty()

    //
    // MODULES: RUN SAMTOOLS VIEW to downsample BAM files
    //
    SAMTOOLS_DOWNSAMPLE (
                ch_bam,
                ch_fasta,
                [], // qname
                'csi' // index_format, for large genomes csi is required
                )
    ch_versions = ch_versions.mix(SAMTOOLS_DOWNSAMPLE.out.versions.first())

    HAPLOTAG (
                SAMTOOLS_DOWNSAMPLE.out.bam
                )
    ch_versions = ch_versions.mix(HAPLOTAG.out.versions.first())

    emit:
    bam      = HAPLOTAG.out.bam      // channel: [ val(meta), path(bam) ]


    versions = ch_versions                    // channel: [ path(versions.yml) ]
}
