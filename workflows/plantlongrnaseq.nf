/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC                                     } from '../modules/nf-core/fastqc/main'
include { MULTIQC                                    } from '../modules/nf-core/multiqc/main'
include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_GENOME     } from '../modules/nf-core/minimap2/align/main'
include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_TRANSCRIPT } from '../modules/nf-core/minimap2/align/main'
include { GFFREAD                                    } from '../modules/nf-core/gffread/main'
include { GUNZIP as GUNZIP_FASTA                     } from '../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_GFF                       } from '../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_GTF                       } from '../modules/nf-core/gunzip/main'
include { GFFREAD as GFFREAD_TRANSCRIPT              } from '../modules/nf-core/gffread/main'
include { paramsSummaryMap                           } from 'plugin/nf-schema'
include { paramsSummaryMultiqc                       } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML                     } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText                     } from '../subworkflows/local/utils_nfcore_plantlongrnaseq_pipeline'
include { BAM_STATS_SAMTOOLS                         } from '../subworkflows/nf-core/bam_stats_samtools/main'
include { CENTRIFUGE_KREPORT                         } from '../modules/nf-core/centrifuge/kreport/main'
include { SAMTOOLS_BAM2FQ                            } from '../modules/nf-core/samtools/bam2fq/main'
include { CENTRIFUGE_CENTRIFUGE                      } from '../modules/nf-core/centrifuge/centrifuge/main'
include { PICARD_FILTERSAMREADS                      } from '../modules/nf-core/picard/filtersamreads/main'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_GENOME      } from '../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_TRANSRIPTOME } from '../modules/nf-core/samtools/sort/main'
include { PBTK_BAM2FASTQ                             } from '../modules/nf-core/pbtk/bam2fastq/main'
include { OARFISH                                    } from '../../modules/modules/nf-core/oarfish/main'
include { DESEQ2_QC                                   } from '../modules/local/deseq2_qc/main'
include { MERGE_COUNTS                                } from '../modules/local/merge_counts'
include { RUN_SQANTI_READS                        } from '../subworkflows/local/run_sqanti_reads'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


// Header files for MultiQC
ch_pca_header_multiqc           = file("$projectDir/workflows/assets/multiqc/deseq2_pca_header.txt", checkIfExists: true)
ch_clustering_header_multiqc    = file("$projectDir/workflows/assets/multiqc/deseq2_clustering_header.txt", checkIfExists: true)

workflow PLANTLONGRNASEQ {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    reads       // channel: path(genome.fasta)
    gff         // channel: path(genome.gff)
    gtf         // channel: path(genome.gtf)

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()


    //
    // Uncompress genome fasta file if required
    //
    if (reads.endsWith('f*q.gz')) {
        ch_fasta    = GUNZIP_FASTA ( [ [:], file(reads, checkIfExists: true) ] ).gunzip.map { it[1] }
        ch_versions = ch_versions.mix(GUNZIP_FASTA.out.versions)
    } else if (reads.endsWith('.bam')) {
        ch_fasta = PBTK_BAM2FASTQ ( [ [:], file(reads, checkIfExists: true) ] ).fasta.map { it[1] }
        ch_versions = ch_versions.mix(PBTK_BAM2FASTQ.out.versions)
    } else if (reads.endsWith('.fasta') || reads.endsWith('.fa')) {
        ch_fasta    = Channel.value(file(reads, checkIfExists: true))
    }


    //
    // Uncompress GTF annotation file or create from GFF3 if required
    //
    if (gtf || gff) {
        if (gtf) {
            if (gtf.endsWith('.gz')) {
                ch_gtf      = GUNZIP_GTF ( [ [:], file(gtf, checkIfExists: true) ] ).gunzip.map { it[1] }
                ch_versions = ch_versions.mix(GUNZIP_GTF.out.versions)
            } else {
                ch_gtf = Channel.value(file(gtf, checkIfExists: true))
            }
        } else if (gff) {
            if (gff.endsWith('.gz')) {
                ch_gff      = GUNZIP_GFF ( [ [:], file(gff, checkIfExists: true) ] ).gunzip
                ch_versions = ch_versions.mix(GUNZIP_GFF.out.versions)
            } else {
                ch_gff = Channel.value(file(gff, checkIfExists: true)).map { [ [:], it ] }
            }
            ch_gtf      = GFFREAD ( ch_gff, [] ).gtf.map { it[1] }
            ch_versions = ch_versions.mix(GFFREAD.out.versions)
        }
    }
    //
    // MODULE: Run FastQC
    //
    FASTQC (
        ch_samplesheet
    )

    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())
    //
    // MODULE: Run GFFREAD to extract transcript sequences
    //
    // extract the spliced transcripts
    GFFREAD_TRANSCRIPT(ch_gtf.map { gtf -> [["id": gtf.simpleName], gtf] },
                       ch_fasta)

    ch_versions = ch_versions.mix(GFFREAD_TRANSCRIPT.out.versions.first())
    //
    // MODULE: Run Minimap2 alignment on gennome
    //
    MINIMAP2_ALIGN_GENOME (
        ch_samplesheet,
        ch_fasta.map { [ [:], it ] },
        true,
        'bai',
        false,
        true
    )

    ch_versions = ch_versions.mix(MINIMAP2_ALIGN_GENOME.out.versions.first())

    //
    // MODULE: Run Minimap2 alignment on transcripts
    //
    MINIMAP2_ALIGN_TRANSCRIPT (
        ch_samplesheet,
        GFFREAD_TRANSCRIPT.out.gffread_fasta,
        true,
        'bai',
        false,
        true
    )

    ch_versions = ch_versions.mix(MINIMAP2_ALIGN_TRANSCRIPT.out.versions.first())
    //
    //MODULE: Run SAMtools sort for transcriptome
    //
    SAMTOOLS_SORT_TRANSRIPTOME (
        MINIMAP2_ALIGN_TRANSCRIPT.out.bam,
        ch_fasta.map { [ [:], it ]}
    )


    //
    // MODULE: Run Oarfish to quantify transcripts aligned
    //
    OARFISH (
        SAMTOOLS_SORT_TRANSRIPTOME.out.bam)

    ch_versions = ch_versions.mix(OARFISH.out.versions.first())



    //
    // MODULE: Merge counts from Oarfish
    //

    combined_ch = OARFISH.out.results
    .toSortedList { a, b ->
        // Sort by sample ID
        a[0].id <=> b[0].id
    }
    .map { tuples ->
        def metas = tuples.collect { it[0] }  // Extract all meta maps
        def files = tuples.collect { it[1] }  // Extract all quant files
        [metas, files]
    }


    MERGE_COUNTS(combined_ch)


    //
    // MODULE: Run DESeq2 QC
    //
    if (!params.skip_deseq2_qc){
        DESEQ2_QC (
            MERGE_COUNTS.out.merged_counts,
            ch_pca_header_multiqc,
            ch_clustering_header_multiqc
        )
        ch_multiqc_files = ch_multiqc_files.mix(DESEQ2_QC.out.pca_multiqc.collect())
        ch_multiqc_files = ch_multiqc_files.mix(DESEQ2_QC.out.dists_multiqc.collect())
        ch_versions = ch_versions.mix(DESEQ2_QC.out.versions.first())

    }




    //
    // SUBWORKFLOW: Run SAMtools stats, flagstat and idxstats
    //

    BAM_STATS_SAMTOOLS ( MINIMAP2_ALIGN_GENOME.out.bam.join(MINIMAP2_ALIGN_GENOME.out.index),
                         ch_fasta.map { [ [:], it ] }
                        )

    ch_multiqc_files = ch_multiqc_files.mix(BAM_STATS_SAMTOOLS.out.stats.collect{it[1]})
    ch_multiqc_files = ch_multiqc_files.mix(BAM_STATS_SAMTOOLS.out.flagstat.collect{it[1]})
    ch_multiqc_files = ch_multiqc_files.mix(BAM_STATS_SAMTOOLS.out.idxstats.collect{it[1]})
    ch_versions = ch_versions.mix(BAM_STATS_SAMTOOLS.out.versions.first())

    //
    // MODULE: Run SAMtools sort for picard filtersamreads
    //
    SAMTOOLS_SORT_GENOME (
        MINIMAP2_ALIGN_GENOME.out.bam,
        ch_fasta.map { [ [:], it ]}
    )

    ch_versions = ch_versions.mix(SAMTOOLS_SORT_GENOME.out.versions.first())

    //
    // MODULE: Run filtersamreads to get unaligned reads
    //
    PICARD_FILTERSAMREADS (
        SAMTOOLS_SORT_GENOME.out.bam.join(MINIMAP2_ALIGN_GENOME.out.index),
        "excludeAligned"
    )

    ch_versions = ch_versions.mix(PICARD_FILTERSAMREADS.out.versions.first())

    // MODULE: Run SAMtools bam2fq to convert unaligned reads to fastq
    SAMTOOLS_BAM2FQ (
        PICARD_FILTERSAMREADS.out.bam,
        false
    )

    //
    // MODULE: Run Centrifuge
    //
    CENTRIFUGE_CENTRIFUGE (
        SAMTOOLS_BAM2FQ.out.reads,
        params.centrifuge_db,
        false,
        false
    )

    ch_versions = ch_versions.mix(CENTRIFUGE_CENTRIFUGE.out.versions.first())

    CENTRIFUGE_KREPORT (
        CENTRIFUGE_CENTRIFUGE.out.report,
        params.centrifuge_db
    )

    ch_multiqc_files = ch_multiqc_files.mix(CENTRIFUGE_KREPORT.out.kreport.collect{it[1]})

    //
    // SUBWORKFLOW: Run SQANTI on reads
    //
    RUN_SQANTI_READS(
                    MINIMAP2_ALIGN_GENOME.out.bam.join(MINIMAP2_ALIGN_GENOME.out.index),
                    ch_fasta.map { [ [:], it ] },
                    ch_gtf
                    )

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  + 'pipeline_software_' +  'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
