### Atlantic
nextflow  run main.nf -resume -profile singularity --input assets/samplesheet_atlantic.csv --outdir output_atlantic_cds --fasta /scratch/nadjafn/reference/Atlantic/ATL_v3.asm.fa  --gtf /scratch/nadjafn/reference/Atlantic/ATL_v3.hc_gene_models.repr.longest_isoforms_exon2cds150.gtf   --centrifuge_db /biodbs/centrifuge/dbs_v2018/ --sqanti_dir /scratch/nadjafn/sqanti3/release_sqanti3 -bg



 ##### Desiree with original phased deiree annotation longest isoform no UTRs
nextflow run main.nf -resume -profile singularity \
                    --input assets/samplesheet_desiree.csv \
                    --outdir output_desiree \
                    --fasta /scratch/nadjafn/reference/Desiree_v1/De_v1_no_scaffold_chloroplast_mt.fa \
                    --gtf /scratch/nadjafn/reference/Desiree_v1/De_v1.functional_annotations_nodigits.longest_isoforms_cds2exon150.gtf \
                    --centrifuge_db /biodbs/centrifuge/dbs_v2018/ \
                    --sqanti_dir /scratch/nadjafn/sqanti3/release_sqanti3 \
                    --sqanti_test true


#### Desiree reads with Unitato reference
nextflow run main.nf -resume -profile singularity \
                    --input assets/samplesheet_desiree.csv \
                    --outdir output_desiree_unitato \
                    --fasta /scratch/markop/WORK/_p_Single-cell/_I_Optimization_of_protocols/_S_Chromium/_A_03_CellRanger-dry/input/UniTato_nuc_mt_plast.fa  \
                    --gtf  /scratch/nadjafn/reference/ADAPT_liftoff/UniTato_nuc_mt_plast.manualfix.nfn.gtf \
                    --centrifuge_db /biodbs/centrifuge/dbs_v2018/ \
                    --sqanti_dir /scratch/nadjafn/sqanti3/release_sqanti3 \
                    --sqanti_test -bg



### Desiree with desiree liftoff reference
nextflow run main.nf -resume -profile singularity \
                    --input assets/samplesheet_desiree.csv \
                    --outdir output_desiree_liftoff_phased \
                    --fasta /scratch/nadjafn/reference/Desiree_v1/De_v1_no_scaffold_chloroplast_mt.fa \
                    --gtf /scratch/nadjafn/reference/Desiree_v1/De_v1.unitato_liftoff_haplotap_gffread.with_chloroplast_and_mito.gtf \
                    --centrifuge_db /biodbs/centrifuge/dbs_v2018/ \
                    --sqanti_dir /scratch/nadjafn/SQANTI3-5.1.2/SQANTI3 \
                    -bg


# add the chloroplast and mitochondria to the annotation

### Desiree

### gff3 to gtf, use agat to perseve genes for bambu
agat_convert_sp_gff2gtf.pl --gff /scratch/nadjafn/reference/Desiree_v1/De_v1.unitato_liftoff_haplotap.gff3 -o /scratch/nadjafn/reference/De_v1.unitato_liftoff_haplotap.agat.gtf

grep "OR9" /scratch/markop/WORK/_p_Single-cell/_I_Optimization_of_protocols/_S_Chromium/_A_03_CellRanger-dry/input/UniTato_nuc_mt_plast_PVY.manualfix.gtf > /scratch/nadjafn/reference/organelles/potato_mito.gtf

grep "NC_" /scratch/markop/WORK/_p_Single-cell/_I_Optimization_of_protocols/_S_Chromium/_A_03_CellRanger-dry/input/UniTato_nuc_mt_plast_PVY.manualfix.gtf > /scratch/nadjafn/reference/organelles/potato_chloroplast.gtf

cat   /scratch/nadjafn/reference/De_v1.unitato_liftoff_haplotap.agat.gtf /scratch/nadjafn/reference/organelles/potato_chloroplast.gtf /scratch/nadjafn/reference/organelles/potato_mito.gtf > /scratch/nadjafn/reference/Desiree_v1/De_v1.unitato_liftoff_haplotap_gffread.with_chloroplast_and_mito.gtf


# run Atlantic with all isoforms but no UTRs
nextflow run main.nf -resume -profile singularity \
                    --input assets/samplesheet_atlantic.csv \
                    --outdir output_atlantic_all_isoforms_no_UTR \
                    --fasta /scratch/nadjafn/reference/Atlantic/ATL_v3.asm.with_chloroplast_and_mito.fa \
                    --gtf  /scratch/nadjafn/reference/Atlantic/ATL_v3.hc_gene_models.agat.exon2cds150.with_chloroplast_and_mito.no_scaffold.gtf \
                    --centrifuge_db /biodbs/centrifuge/dbs_v2018/ \
                    --sqanti_dir /scratch/nadjafn/sqanti3/release_sqanti3 \
                    --sqanti_test -bg --technology ONT
