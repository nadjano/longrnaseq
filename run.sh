 nextflow  run main.nf -resume -profile singularity --input assets/samplesheet_atlantic.csv --outdir output --fasta /scratch/nadjafn/reference/Atlantic/ATL_v3.asm.fa  --gtf /scratch/nadjafn/reference/Atlantic/ATL_v3.hc_gene_models.repr.gtf  --centrifuge_db /biodbs/centrifuge/dbs_v2018/ --skip_deseq2_qc --sqanti_dir /scratch/nadjafn/sqanti3/release_sqanti3 -bg



 ##### Desiree
nextflow run main.nf -resume -profile singularity \
                    --input assets/samplesheet_desiree.csv \
                    --outdir output_desiree \
                    --fasta /scratch/nadjafn/reference/Desiree_v1/De_v1_no_scaffold_chloroplast_mt.fa \
                    --gtf /scratch/nadjafn/reference/Desiree_v1/De_v1.functional_annotations_nodigits.longest_isoforms_cds2exon150.gtf \
                    --centrifuge_db /biodbs/centrifuge/dbs_v2018/ \
                    --sqanti_dir /scratch/nadjafn/sqanti3/release_sqanti3 \
                    --sqanti_test true


#### Desiree with Unitato reference
nextflow run main.nf -resume -profile singularity \
                    --input assets/samplesheet_desiree.csv \
                    --outdir output_desiree_unitato \
                    --fasta /scratch/markop/WORK/_p_Single-cell/_I_Optimization_of_protocols/_S_Chromium/_A_03_CellRanger-dry/input/UniTato_nuc_mt_plast.fa  \
                    --gtf /scratch/markop/WORK/_p_Single-cell/_I_Optimization_of_protocols/_S_Chromium/_A_03_CellRanger-dry/input/UniTato_nuc_mt_plast.manualfix.gtf  \
                    --centrifuge_db /biodbs/centrifuge/dbs_v2018/ \
                    --sqanti_dir /scratch/nadjafn/sqanti3/release_sqanti3 \
                    --sqanti_test -bg
