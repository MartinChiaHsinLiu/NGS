
conda activate IO
cd /mnt/ngs_data/liu/NDMC/SC0205/
#GATK 可利用bioconda來安裝
#samtools 可利用bioconda來安裝 conda install -c bioconda samtools
#vcftools 可利用bioconda來安裝 conda install -c bioconda vcftools
#hisat2 可利用bioconda來安裝 conda install hisat2
#stringtie 可利用bioconda來安裝 conda install -c bioconda stringtie
export PERL5LIB=/mnt/ngs_data/liu/ensembl-vep
export PERL5LIB=/mnt/ngs_data/liu/vcftools_0.1.13/perl/

#flexbar用sudo apt install flexbar安裝
#flexbar -n 16 -r Tumor_RNA_R1.fastq.gz -p Tumor_RNA_R2.fastq.gz --adapter-trim-end RIGHT --zip-output GZ --adapter-preset TruSeq -ap ON #do not work in openstack due to no preset adapter
flexbar -n 16 -r Tumor_RNA_R1.fastq.gz -p Tumor_RNA_R2.fastq.gz --adapter-trim-end RIGHT --zip-output GZ -a /mnt/ngs_data/liu/grch37_tran/IlluminaTruSeqAdaptor.fasta 

hisat2 -p 16 -x /mnt/ngs_data/liu/grch37_tran/genome_tran -1 flexbarOut_1.fastq.gz -2 flexbarOut_2.fastq.gz -S Tumor_RNA.sam
samtools view -bS Tumor_RNA.sam > Tumor_RNA.bam
samtools sort -o Tumor_RNA_sorted.bam Tumor_RNA.bam
samtools index Tumor_RNA_sorted.bam

gatk AddOrReplaceReadGroups -I Tumor_RNA_sorted.bam -O rg_added_sorted.bam -SO coordinate -RGID Tumor_RNA -RGLB TruSeq -RGPL illumina -RGPU novaseq -RGSM Tumor_RNA
gatk MarkDuplicates -I rg_added_sorted.bam -O dedupped.bam  -CREATE_INDEX true -VALIDATION_STRINGENCY SILENT -M output.metrics 
gatk BaseRecalibrator -R /mnt/ngs_data/liu/grch37_tran/Homo_sapiens.GRCh37.75.dna.chromosome.all.fasta -I dedupped.bam -O Tumor_RNA.recal_data_before.grp --known-sites /home/yaron/appreci8/resources/dbsnp_138.b37.vcf --known-sites /home/yaron/appreci8/resources/Mills_and_1000G_gold_standard.indels.b37.vcf --known-sites /home/yaron/appreci8/resources/1000G_phase1.indels.b37.vcf 
gatk ApplyBQSR -R /mnt/ngs_data/liu/grch37_tran/Homo_sapiens.GRCh37.75.dna.chromosome.all.fasta -I dedupped.bam -O Tumor_RNA.recali.bam -bqsr Tumor_RNA.recal_data_before.grp
gatk SplitNCigarReads -R /mnt/ngs_data/liu/grch37_tran/Homo_sapiens.GRCh37.75.dna.chromosome.all.fasta -I Tumor_RNA.recali.bam -O split.bam
gatk HaplotypeCaller -R /mnt/ngs_data/liu/grch37_tran/Homo_sapiens.GRCh37.75.dna.chromosome.all.fasta -I split.bam --dont-use-soft-clipped-bases --dbsnp /home/yaron/appreci8/resources/dbsnp_138.b37.vcf -stand-call-conf 20 -O Tumor_RNA_sorted.HC.vcf
gatk VariantFiltration -R /mnt/ngs_data/liu/grch37_tran/Homo_sapiens.GRCh37.75.dna.chromosome.all.fasta -V Tumor_RNA_sorted.HC.vcf -window 35 -cluster 3 --filter-name FS --filter-expression "FS > 30.0" --filter-name QD --filter-expression "QD < 2.0" -O Tumor_RNA_sorted.raw.snps.indels.vcf 

#/mnt/ngs_data/liu/grch37_tran/Homo_sapiens.GRCh37.75_autosome.gff3 需要抓取, 這是去掉mitochondria的 gene file format version 3 reference
stringtie Tumor_RNA_sorted.bam -G /mnt/ngs_data/liu/grch37_tran/Homo_sapiens.GRCh37.75_autosome.gff3 -o Tumor_RNA_sorted.stringtie.gff -A gene_abund.out -C transcript_abund.gtf

python /mnt/ngs_data/liu/gtf_to_csv.py -i Tumor_RNA_sorted.stringtie.gtf -o Tumor_RNA_sorted.stringtie.csv

awk -F, 'gsub(/,/,"\t")1' Tumor_RNA_sorted.stringtie.csv > Tumor_RNA_sorted.stringtie.tsv






