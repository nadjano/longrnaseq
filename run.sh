 nextflow  run main.nf -resume -profile singularity --input assets/samplesheet_atlantic.csv --outdir output --fasta /scratch/nadjafn/reference/Atlantic/ATL_v3.asm.fa  --gtf /scratch/nadjafn/reference/Atlantic/ATL_v3.hc_gene_models.repr.gtf  --centrifuge_db /biodbs/centrifuge/dbs_v2018/



 ##### Desiree
nextflow run main.nf -resume -profile singularity \
                    --input assets/samplesheet_desiree.csv \
                    --outdir output_desiree \
                    --fasta /scratch/nadjafn/reference/Desiree_v1/De_v1_no_scaffold_chloroplast_mt.fa \
                    --gtf /scratch/nadjafn/reference/Desiree_v1/De_v1.functional_annotations_nodigits.longest_isoforms_cds150.gtf \
                    --centrifuge_db /biodbs/centrifuge/dbs_v2018/ \
                    -bg