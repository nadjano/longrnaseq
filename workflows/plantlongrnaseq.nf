/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC                                     } from '../modules/nf-core/fastqc/main'
include { MULTIQC                                    } from '../modules/nf-core/multiqc/main'
include { MINIMAP2_ALIGN                             } from '../modules/nf-core/minimap2/align/main'
include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_TRANSCRIPT } from '../modules/nf-core/minimap2/align/main'
include { GFFREAD                                    } from '../modules/nf-core/gffread/main'
include { GFFREAD as GFFREAD_TRANSCRIPT              } from '../modules/nf-core/gffread/main'
include { paramsSummaryMap                           } from 'plugin/nf-schema'
include { paramsSummaryMultiqc                       } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML                     } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText                     } from '../subworkflows/local/utils_nfcore_plantlongrnaseq_pipeline'
include { BAM_STATS_SAMTOOLS                         } from '../subworkflows/nf-core/bam_stats_samtools/main'   
include { KRAKEN2_KRAKEN2                            } from '../modules/nf-core/kraken2/kraken2/main'  
include { PICARD_FILTERSAMREADS                      } from '../modules/nf-core/picard/filtersamreads/main' 

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow PLANTLONGRNASEQ {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    fasta       // channel: path(genome.fasta)
    gff         // channel: path(genome.gff)
    gtf         // channel: path(genome.gtf)
    
    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()


     //
    // Uncompress genome fasta file if required
    //
    if (fasta.endsWith('.gz')) {
        ch_fasta    = GUNZIP_FASTA ( [ [:], file(fasta, checkIfExists: true) ] ).gunzip.map { it[1] }
        ch_versions = ch_versions.mix(GUNZIP_FASTA.out.versions)
    } else {
        ch_fasta = Channel.value(file(fasta, checkIfExists: true))
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

    GFFREAD_TRANSCRIPT.out.gffread_fasta.view()

    ch_versions = ch_versions.mix(GFFREAD_TRANSCRIPT.out.versions.first())
    //
    // MODULE: Run Minimap2 alignment on gennome
    //
    MINIMAP2_ALIGN (
        ch_samplesheet,
        ch_fasta.map { [ [:], it ] },
        true,
        'bai',
        false,
        true
    )
    MINIMAP2_ALIGN.out.bam.view()
    MINIMAP2_ALIGN.out.index.view()

    ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions.first())
    
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
    // SUBWORKFLOW: Run SAMtools stats, flagstat and idxstats
    //

    BAM_STATS_SAMTOOLS ( MINIMAP2_ALIGN.out.bam.join(MINIMAP2_ALIGN.out.index), 
                         ch_fasta.map { [ [:], it ] }
                        )

    ch_multiqc_files = ch_multiqc_files.mix(BAM_STATS_SAMTOOLS.out.stats.collect{it[1]})
    ch_versions = ch_versions.mix(BAM_STATS_SAMTOOLS.out.versions.first())

    //
    // MODULE: Run filtersamreads to get unaligned reads
    //
    PICARD_FILTERSAMREADS (
        MINIMAP2_ALIGN.out.bam.join(MINIMAP2_ALIGN.out.index),
        "excludeAligned"
    )


    // 
    // MODULE: Run Kraken2
    //
//     KRAKEN2_KRAKEN (
//             ch_unaligned_sequences,
//             params.kraken_db,
//             params.save_kraken_assignments,
//             params.save_kraken_unassigned        
// )

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
