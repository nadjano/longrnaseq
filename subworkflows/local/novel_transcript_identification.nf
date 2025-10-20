// 
// Subworkflow for novel transcript identification using BAMBU, LIFTOFF, and GFFCOMPARE
//

include { BAMBU                                        } from '../../modules/local/bambu'
include { GAWK                                         } from '../../modules/nf-core/gawk/main'
include { SQANTIQC                                     } from '../../modules/local/sqanti/sqantiqc'
include { SAMTOOLS_VIEW as SAMTOOLS_FILTER             } from '../../modules/nf-core/samtools/view/main'
include { LIFTOFF                                      } from '../../modules/nf-core/liftoff/main'
include { GFFCOMPARE                                   } from '../../modules/nf-core/gffcompare/main'
include { GFF_ATTR_RENAMING                            } from '../../modules/local/gff_attr_renaming/main'
workflow NOVEL_TRANSCRIPT_IDENTIFICATION {

    take:
    ch_bam          // channel: [ val(meta), path(bam) ]
    ch_fasta        // channel: [ path(fasta) ]
    ch_gtf          // channel:  path(gtf)
    
    main:

    ch_versions = Channel.empty()

    //
    // Prepare channels for BAMBU
    //
    
    // Group BAM files by meta for BAMBU (which expects multiple BAM files)
    ch_bam_grouped = ch_bam.toSortedList { a, b -> a[0].id <=> b[0].id }.map { it.collect { tuple -> tuple[1] } }
    // Prepare GTF with simplified meta for BAMBU
    ch_gtf_for_bambu = ch_gtf.map { gtf -> [["id": gtf.simpleName], gtf] }

    //
    // MODULE: Run BAMBU to reconstruct transcripts
    //
    ch_fasta.view()
    BAMBU(
        ch_fasta,
        ch_gtf_for_bambu,
        ch_bam_grouped
    )
    ch_versions = ch_versions.mix(BAMBU.out.versions.first())

    // 
    // MODULE: Filter BAMBU annotation with GAWK
    //
    GAWK(
        BAMBU.out.extended_gtf.map { gtf -> [["id": gtf.simpleName], gtf] },
        [],  // No additional input files
        []   // No program file
    )
    ch_versions = ch_versions.mix(GAWK.out.versions.first())

    //
    // MODULE: Liftoff novel BAMBU isoforms to other haplotypes
    //
    LIFTOFF (
        ch_fasta.map { fasta -> [["id": fasta.simpleName], fasta]  },
        ch_fasta, // just get the reference fasta from id, fasta channel
        GAWK.out.output.map { it[1] },
        []    
    )
    ch_versions = ch_versions.mix(LIFTOFF.out.versions.first())


    // Create a channel with placeholder values
    ch_fasta_empty = Channel.of([ [:], [], [] ])  // Empty meta, empty file list


    ch_combined_gtf = LIFTOFF.out.gff3
    .map { meta, gtf -> gtf.copyTo("liftoff.gtf") }
    .combine(
        BAMBU.out.extended_gtf.map { gtf -> gtf.copyTo("bambu.gtf") }
    )
    .map { gtf1, gtf2 -> [[id: "combined"], [gtf1, gtf2]] }

    ch_combined_gtf.view()
    // MODDULE: Merge liftoff with previous annotation
    GFFCOMPARE ( 
        ch_combined_gtf,
        ch_fasta_empty,  // tuple val(meta2), path(fasta), path(fai)
        BAMBU.out.extended_gtf.map { gtf -> [["id": gtf.simpleName], gtf] }

    )
    ch_versions = ch_versions.mix(GFFCOMPARE.out.versions.first())


    //
    // MODULE: Fix gtf file gene_ids and transcript_ids as gffcompare changes them
    //
    GFF_ATTR_RENAMING (
        GFFCOMPARE.out.combined_gtf
    )
    ch_versions = ch_versions.mix(GFF_ATTR_RENAMING.out.versions.first())

    emit:
    // TODO nf-core: edit emitted channels

    versions = ch_versions
    novel_gtf =  GFF_ATTR_RENAMING.out.fixed_gtf


}


