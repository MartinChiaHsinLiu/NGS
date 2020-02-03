
conda activate IO
cd /mnt/ngs_data/liu/NDMC/SC0205/
#GATK4 可利用bioconda來安裝 conda install -c bioconda gatk4
#ensembl-vep 可利用bioconda來安裝 conda install -c bioconda ensembl-vep
#vcftools 可利用bioconda來安裝 conda install -c bioconda vcftools
#samtools 可利用bioconda來安裝 conda install -c bioconda samtools
#bcftools 可利用bioconda來安裝 conda install -c bioconda bcftools
#vt 可利用bioconda來安裝 conda install -c bioconda vt
#vatools 可利用pip來安裝 pip install vatools
#pvactools 可利用pip來安裝 pip install pvactools
export PERL5LIB=/mnt/ngs_data/liu/ensembl-vep
export PERL5LIB=/mnt/ngs_data/liu/vcftools_0.1.13/perl/


gatk LeftAlignAndTrimVariants -R /mnt/ngs_data/liu/grch37_tran/Homo_sapiens.GRCh37.75.dna.chromosome.all.fasta --variant SC0205_test.recali-filtered.vcf -O SC0205_test.mutect2_LAATV.vcf
bgzip SC0205_test.mutect2_LAATV.vcf
tabix -p vcf SC0205_test.mutect2_LAATV.vcf.gz

vt decompose -s SC0205_test.mutect2_LAATV.vcf.gz -o SC0205_test.recali-filtered.decomposed.vcf
gatk IndexFeatureFile -F SC0205_test.recali-filtered.decomposed.vcf
bgzip SC0205_test.recali-filtered.decomposed.vcf
tabix -p vcf SC0205_test.recali-filtered.decomposed.vcf.gz

vt decompose -s Tumor_RNA_sorted.raw.snps.indels.vcf -o Tumor_RNA_sorted.raw.snps.indels.decomposed.vcf
gatk IndexFeatureFile -F Tumor_RNA_sorted.raw.snps.indels.decomposed.vcf
bgzip Tumor_RNA_sorted.raw.snps.indels.decomposed.vcf
tabix -p vcf Tumor_RNA_sorted.raw.snps.indels.decomposed.vcf.gz



bcftools isec SC0205_test.recali-filtered.decomposed.vcf.gz Tumor_RNA_sorted.raw.snps.indels.decomposed.vcf.gz -p intersect
cd intersect
bcftools annotate -x ^FORMAT/GT,FORMAT/AD,FORMAT/DP 0002.vcf > 0002_trim.vcf
bcftools annotate -x ^FORMAT/GT,FORMAT/AD,FORMAT/DP 0003.vcf > 0003_trim.vcf
bgzip 0002_trim.vcf
bgzip 0003_trim.vcf
tabix -p vcf 0002_trim.vcf.gz
tabix -p vcf 0003_trim.vcf.gz
bcftools merge 0002_trim.vcf.gz 0003_trim.vcf.gz > merge.vcf

cp merge.vcf /mnt/ngs_data/liu/NDMC/SC0205/
cd ..
cat merge.vcf | vcf-sort > SC0205_test.DNARNA.recali-filtered.vcf
bgzip SC0205_test.DNARNA.recali-filtered.vcf
tabix -p vcf SC0205_test.DNARNA.recali-filtered.vcf.gz

#這3步可以平行
python /mnt/ngs_data/liu/bin/bam_readcount_helper.py SC0205_test.DNARNA.recali-filtered.vcf.gz Tumor_RNA /mnt/ngs_data/liu/grch37_tran/Homo_sapiens.GRCh37.75.dna.chromosome.all.fasta Tumor_RNA.recali.bam /mnt/ngs_data/liu/NDMC/SC0205/coverage/
python /mnt/ngs_data/liu/bin/bam_readcount_helper.py SC0205_test.DNARNA.recali-filtered.vcf.gz Tumor_WES /mnt/ngs_data/liu/grch37_tran/Homo_sapiens.GRCh37.75.dna.chromosome.all.fasta Tumor_WES.recali.bam /mnt/ngs_data/liu/NDMC/SC0205/coverage/
python /mnt/ngs_data/liu/bin/bam_readcount_helper.py SC0205_test.DNARNA.recali-filtered.vcf.gz Normal_WES /mnt/ngs_data/liu/grch37_tran/Homo_sapiens.GRCh37.75.dna.chromosome.all.fasta Normal_WES.recali.bam /mnt/ngs_data/liu/NDMC/SC0205/coverage/




vcf-readcount-annotator SC0205_test.DNARNA.recali-filtered.vcf.gz /mnt/ngs_data/liu/NDMC/SC0205/coverage/Tumor_RNA_bam_readcount_snv.tsv RNA -s Tumor_WES -t snv -o SC0205_test.DNARNA.recali-filtered.RNAsnv.vcf

vcf-readcount-annotator SC0205_test.DNARNA.recali-filtered.RNAsnv.vcf /mnt/ngs_data/liu/NDMC/SC0205/coverage/Tumor_RNA_bam_readcount_indel.tsv RNA -s Tumor_WES -t indel -o SC0205_test.DNARNA.recali-filtered.RNAsnv.RNAindel.vcf

vcf-readcount-annotator SC0205_test.DNARNA.recali-filtered.RNAsnv.RNAindel.vcf /mnt/ngs_data/liu/NDMC/SC0205/coverage/Tumor_WES_bam_readcount_snv.tsv DNA -s Tumor_WES -t snv -o SC0205_test.DNARNA.recali-filtered.RNAsnv.RNAindel.Tsnv.vcf

vcf-readcount-annotator SC0205_test.DNARNA.recali-filtered.RNAsnv.RNAindel.Tsnv.vcf /mnt/ngs_data/liu/NDMC/SC0205/coverage/Tumor_WES_bam_readcount_indel.tsv DNA -s Tumor_WES -t indel -o SC0205_test.DNARNA.recali-filtered.RNAsnv.RNAindel.Tsnv.Tindel.vcf

vcf-readcount-annotator SC0205_test.DNARNA.recali-filtered.RNAsnv.RNAindel.Tsnv.Tindel.vcf /mnt/ngs_data/liu/NDMC/SC0205/coverage/Normal_WES_bam_readcount_snv.tsv DNA -s Normal_WES -t snv -o SC0205_test.DNARNA.recali-filtered.RNAsnv.RNAindel.Tsnv.Tindel.Nsnv.vcf

vcf-readcount-annotator SC0205_test.DNARNA.recali-filtered.RNAsnv.RNAindel.Tsnv.Tindel.Nsnv.vcf /mnt/ngs_data/liu/NDMC/SC0205/coverage/Normal_WES_bam_readcount_indel.tsv DNA -s Normal_WES -t indel -o SC0205_test.DNARNA.recali-filtered.RNAsnv.RNAindel.Tsnv.Tindel.Nsnv.Nindel.vcf


vep --input_file SC0205_test.DNARNA.recali-filtered.RNAsnv.RNAindel.Tsnv.Tindel.Nsnv.Nindel.vcf --output_file SC0205_test.DNARNA.recali-filtered.RNAsnv.RNAindel.Tsnv.Tindel.Nsnv.Nindel.vep.vcf --format vcf --vcf --symbol --terms SO --tsl --hgvs --fasta /mnt/ngs_data/liu/grch37_tran/Homo_sapiens.GRCh37.75.dna.chromosome.all.fasta --offline --cache --dir_cache /mnt/ngs_data/liu/grch37_tran --plugin Downstream --plugin Wildtype --dir_plugins /mnt/ngs_data/liu/grch37_tran/VEP_plugins/ -o SC0205_test.DNARNA.recali-filtered.RNAsnv.RNAindel.Tsnv.Tindel.Nsnv.Nindel.vep.vcf

awk -F '\t' '{if($0 ~ /\#/) print; else if($7 == "PASS") print}' filtered.vcf > PASS.vcf
vep --input_file PASS.vcf --output_file PASS.vep.vcf --format vcf --vcf --symbol --terms SO --tsl --hgvs --fasta /mnt/ngs_data/liu/grch37_tran/Homo_sapiens.GRCh37.75.dna.chromosome.all.fasta --offline --cache --dir_cache /mnt/ngs_data/liu/grch37_tran --plugin Downstream --plugin Wildtype --dir_plugins /mnt/ngs_data/liu/grch37_tran/VEP_plugins/ 

awk -F '\t' '{if($0 ~ /\#/) print; else if($7 == "PASS") print}' SC0205_test.DNARNA.recali-filtered.RNAsnv.RNAindel.Tsnv.Tindel.Nsnv.Nindel.vcf > PASS_selected.vcf
vep --input_file PASS_selected.vcf --output_file PASS_selected.vep.vcf --format vcf --vcf --symbol --terms SO --tsl --hgvs --fasta /mnt/ngs_data/liu/grch37_tran/Homo_sapiens.GRCh37.75.dna.chromosome.all.fasta --offline --cache --dir_cache /mnt/ngs_data/liu/grch37_tran --plugin Downstream --plugin Wildtype --dir_plugins /mnt/ngs_data/liu/grch37_tran/VEP_plugins/ 


vcf-expression-annotator SC0205_test.DNARNA.recali-filtered.RNAsnv.RNAindel.Tsnv.Tindel.Nsnv.Nindel.vep.vcf gene_abund.out -i 'Gene ID' -e 'TPM' -s Tumor_WES stringtie gene -o SC0205_test.DNARNA.recali-filtered.RNAsnv.RNAindel.Tsnv.Tindel.Nsnv.Nindel.vep.gx.vcf

vcf-expression-annotator --ignore-transcript-version SC0205_test.DNARNA.recali-filtered.RNAsnv.RNAindel.Tsnv.Tindel.Nsnv.Nindel.vep.gx.vcf Tumor_RNA_sorted.stringtie.gff stringtie transcript -s Tumor_WES -o SC0205_test.DNARNA.recali-filtered.RNAsnv.RNAindel.Tsnv.Tindel.Nsnv.Nindel.vep.tx.vcf


samtools sort -T SC0205_test.recali.m2.Tumor_WES.bam -o SC0205_test.recali.m2.Tumor_WES_sorted.bam SC0205_test.recali.m2.Tumor_WES.bam
samtools sort -T SC0205_test.recali.m2.Normal_WES.bam -o SC0205_test.recali.m2.Normal_WES_sorted.bam SC0205_test.recali.m2.Normal_WES.bam

samtools index SC0205_test.recali.m2.Tumor_WES_sorted.bam
samtools index SC0205_test.recali.m2.Normal_WES_sorted.bam
samtools view -bS Tumor_WES_sorted.sam > Tumor_WES_sorted.bam
samtools view -bS Normal_WES_sorted.sam > Normal_WES_sorted.bam
samtools index Tumor_WES_sorted.bam
samtools index Normal_WES_sorted.bam

pvacseq run SC0205_test.DNARNA.recali-filtered.RNAsnv.RNAindel.Tsnv.Tindel.Nsnv.Nindel.vep.tx.vcf Tumor_WES HLA-A*24:02,HLA-A*11:01,HLA-B*15:01,HLA-B*58:01,C*03:02 NetMHC /mnt/ngs_data/liu/NDMC/SC0205/output/SNV/ -e 8,9,10,11 --pass-only --normal-sample-name Normal_WES

#這裡需要先跑完gene_fusion_pipeline.txt
pvacfuse run /mnt/ngs_data/liu/NDMC/SC0205/agfusion/ Tumor_WES HLA-A*24:02,HLA-A*11:01,HLA-B*15:01,HLA-B*58:01 NetMHC /mnt/ngs_data/liu/NDMC/SC0205/agfusion_out/ -e 8,9,10,11

#把pvacseq的結果parse成class I model input
python /mnt/ngs_data/cromwell/pvac_to_model/pvacseq_to_model.py /mnt/ngs_data/liu/NDMC/SC0205/output/SNV/MHC_Class_I/Tumor_WES.tsv /mnt/ngs_data/liu/NDMC/SC0205/output/SNV/MHC_Class_I/Tumor_WES.all_epitopes.tsv /mnt/ngs_data/liu/NDMC/SC0205/SC0205_test.DNARNA.recali-filtered.RNAsnv.RNAindel.Tsnv.Tindel.Nsnv.Nindel.vep.tx.vcf /mnt/ngs_data/liu/NDMC/SC0205/output/SNV/MHC_Class_I/Tumor_WES.all_epitopes_flanking2.csv











