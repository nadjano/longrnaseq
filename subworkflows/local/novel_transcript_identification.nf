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
include { SPLIT_HAPLOTYPES                       } from '../../modules/local/split_haplotypes/main'
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
    // MODULE: Split GFF by chromosome to avoid copied transcripts to be discarded in LIFTOFF
    //
    SPLIT_HAPLOTYPES (
        GAWK.out.output  

    )
    SPLIT_HAPLOTYPES.out.split_gffs.view()
    //
    // MODULE: Liftoff novel BAMBU isoforms to other haplotypes
    //
    SPLIT_HAPLOTYPES.out.split_gffs
        .transpose()
        .map { meta, gff ->
            def chr = gff.name.replaceAll(/.*\.bambu\.(.*)\.gff/, '$1')
            [meta + ["id": chr], gff]
        }
        .combine(ch_fasta.map { [it] })  // Add fasta to each
        .map { meta, gff, fasta ->
            [meta, fasta, gff]
        }
        .set { ch_liftoff_input }
    ch_liftoff_input.view()
    LIFTOFF (
        ch_liftoff_input.map { meta, fasta, gff -> [meta, fasta] },
        ch_fasta,
        ch_liftoff_input.map { meta, fasta, gff -> gff },
        []
    )
    
    // LIFTOFF (
    //     ch_fasta.map { fasta -> [["id": fasta.simpleName], fasta]  },
    //     ch_fasta, // just get the reference fasta from id, fasta channel
    //     GAWK.out.output.map { it[1] },
    //     []    
    // )
    ch_versions = ch_versions.mix(LIFTOFF.out.versions.first())

    

    // Create a channel with placeholder values
    ch_fasta_empty = Channel.of([ [:], [], [] ])  // Empty meta, empty file list


    // Collect all Liftoff chromosome outputs into a single list
    ch_combined_gtf = LIFTOFF.out.gff3
        .map { meta, gtf -> gtf }
        .collect() // Collects all files into a single list
        .map { liftoff_gtfs ->
            [[id: "combined"], liftoff_gtfs]
        }
        .combine(
            BAMBU.out.extended_gtf  // Don't map it, use it directly
        )
        .map { meta, liftoff_list, bambu_gtf ->
            [meta, liftoff_list + [bambu_gtf]]  // Wrap bambu_gtf in a list
        }


    ch_combined_gtf.view()

    // Run GFFCOMPARE with all files
    GFFCOMPARE (
        ch_combined_gtf,
        ch_fasta_empty,
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


