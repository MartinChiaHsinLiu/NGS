
conda activate IO
cd /mnt/ngs_data/liu/NDMC/SC0205/

#trimmomatic 可利用bioconda來安裝 conda install -c bioconda trimmomatic
#GATK4 可利用bioconda來安裝 conda install -c bioconda gatk4
#ensembl-vep 可利用bioconda來安裝 conda install -c bioconda ensembl-vep
#vcftools 可利用bioconda來安裝 conda install -c bioconda vcftools
#samtools 可利用bioconda來安裝 conda install -c bioconda samtools
#multiqc 可利用bioconda來安裝 conda install -c bioconda -c conda-forge multiqc
#fastqc 可利用bioconda來安裝 conda install -c bioconda fastqc
export PERL5LIB=/mnt/ngs_data/liu/ensembl-vep
export PERL5LIB=/mnt/ngs_data/liu/vcftools_0.1.13/perl/

#這兩步可以平行
trimmomatic PE -threads 16 -phred33 Tumor_WES_R1.fastq.gz Tumor_WES_R2.fastq.gz Tumor_WES_R1.paired.trim.fastq.gz Tumor_WES_R1.unpaired.trim.fastq.gz Tumor_WES_R2.paired.trim.fastq.gz Tumor_WES_R2.unpaired.trim.fastq.gz ILLUMINACLIP:/opt/miniconda3/share/trimmomatic-0.36-5/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36
trimmomatic PE -threads 16 -phred33 Normal_WES_R1.fastq.gz Normal_WES_R2.fastq.gz Normal_WES_R1.paired.trim.fastq.gz Normal_WES_R1.unpaired.trim.fastq.gz Normal_WES_R2.paired.trim.fastq.gz Normal_WES_R2.unpaired.trim.fastq.gz ILLUMINACLIP:/opt/miniconda3/share/trimmomatic-0.36-5/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36


#這四步可以平行 
#-o 指定output folder 
fastqc -t 6 Normal_WES_R1.paired.trim.fastq.gz -o /mnt/ngs_data/liu/NDMC/SC0205/
fastqc -t 6 Normal_WES_R2.paired.trim.fastq.gz -o /mnt/ngs_data/liu/NDMC/SC0205/
fastqc -t 6 Tumor_WES_R1.paired.trim.fastq.gz  -o /mnt/ngs_data/liu/NDMC/SC0205/
fastqc -t 6 Tumor_WES_R2.paired.trim.fastq.gz  -o /mnt/ngs_data/liu/NDMC/SC0205/

#這兩步可以平行
bwa mem -L 1 -M -t 16 /mnt/ngs_data/liu/grch37_tran/Homo_sapiens.GRCh37.75.dna.chromosome.all.fasta Tumor_WES_R1.paired.trim.fastq.gz Tumor_WES_R2.paired.trim.fastq.gz > Tumor_WES.sam
bwa mem -L 1 -M -t 16 /mnt/ngs_data/liu/grch37_tran/Homo_sapiens.GRCh37.75.dna.chromosome.all.fasta Normal_WES_R1.paired.trim.fastq.gz Normal_WES_R2.paired.trim.fastq.gz > Normal_WES.sam

#這兩步可以平行
samtools sort Tumor_WES.sam -O sam -o Tumor_WES_sorted.sam
samtools sort Normal_WES.sam -O sam -o Normal_WES_sorted.sam

#這兩步可以平行
gatk AddOrReplaceReadGroups -I Tumor_WES_sorted.sam -ID Tumor_WES -LB Tumor_WES -PL illumina -PU Tumor_WES_R1 -SM Tumor_WES -O Tumor_WES_sorted_ARG.bam
gatk AddOrReplaceReadGroups -I Normal_WES_sorted.sam -ID Normal_WES -LB Normal_WES -PL illumina -PU Normal_WES_R1 -SM Normal_WES -O Normal_WES_sorted_ARG.bam
#這兩步可以平行
gatk MarkDuplicates -I Tumor_WES_sorted_ARG.bam -O Tumor_WES.markdup.bam -M Tumor_WES.markdup.metrics
gatk MarkDuplicates -I Normal_WES_sorted_ARG.bam -O Normal_WES.markdup.bam -M Normal_WES.markdup.metrics

#這兩步可以平行
gatk BaseRecalibrator -R /mnt/ngs_data/liu/grch37_tran/Homo_sapiens.GRCh37.75.dna.chromosome.all.fasta -I Tumor_WES.markdup.bam -O Tumor_WES.recal_data_before.grp --known-sites /home/yaron/appreci8/resources/dbsnp_138.b37.vcf --known-sites /home/yaron/appreci8/resources/Mills_and_1000G_gold_standard.indels.b37.vcf --known-sites /home/yaron/appreci8/resources/1000G_phase1.indels.b37.vcf 
gatk BaseRecalibrator -R /mnt/ngs_data/liu/grch37_tran/Homo_sapiens.GRCh37.75.dna.chromosome.all.fasta -I Normal_WES.markdup.bam -O Normal_WES.recal_data_before.grp --known-sites /home/yaron/appreci8/resources/dbsnp_138.b37.vcf --known-sites /home/yaron/appreci8/resources/Mills_and_1000G_gold_standard.indels.b37.vcf --known-sites /home/yaron/appreci8/resources/1000G_phase1.indels.b37.vcf 

#這兩步可以平行
gatk ApplyBQSR -R /mnt/ngs_data/liu/grch37_tran/Homo_sapiens.GRCh37.75.dna.chromosome.all.fasta -I Tumor_WES.markdup.bam -O Tumor_WES.recali.bam -bqsr Tumor_WES.recal_data_before.grp
gatk ApplyBQSR -R /mnt/ngs_data/liu/grch37_tran/Homo_sapiens.GRCh37.75.dna.chromosome.all.fasta -I Normal_WES.markdup.bam -O Normal_WES.recali.bam -bqsr Normal_WES.recal_data_before.grp

#這兩步可以平行
samtools depth Tumor_WES.recali.bam -aa -b /mnt/ngs_data/liu/grch37_tran/hg19.exons.bed > Tumor_WES.depth.txt
samtools depth Normal_WES.recali.bam -aa -b /mnt/ngs_data/liu/grch37_tran/hg19.exons.bed > Normal_WES.depth.txt
#TruSeq在 /mnt/ngs_data/liu/grch37_tran/truseq-exome-targeted-regions-manifest-v1-2.GRCh37.bed
#hg19在 /mnt/ngs_data/liu/grch37_tran/hg19.exons.bed
#Exome_Agilent_V6在 /mnt/ngs_data/liu/grch37_tran/Exome_Agilent_V6.bed

gatk SplitIntervals -R /mnt/ngs_data/liu/grch37_tran/Homo_sapiens.GRCh37.75.dna.chromosome.all.fasta -L /mnt/ngs_data/liu/grch37_tran/Exome_Agilent_V6.bed -scatter 10 -O scattered.intervals

#這兩步可以平行
gatk GetSampleName -I Normal_WES.recali.bam -O normal_name.txt
gatk GetSampleName -I Tumor_WES.recali.bam -O tumor_name.txt


gatk Mutect2 -R /mnt/ngs_data/liu/grch37_tran/Homo_sapiens.GRCh37.75.dna.chromosome.all.fasta -I Normal_WES.recali.bam -tumor Normal_WES --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter -O pon.vcf.gz
gatk Mutect2 --tumor-sample Tumor_WES --normal-sample Normal_WES -I Tumor_WES.recali.bam -I Normal_WES.recali.bam --panel-of-normals pon.vcf.gz --germline-resource /mnt/ngs_data/liu/grch37_tran/af-only-gnomad.raw.sites.b37.vcf.gz --reference /mnt/ngs_data/liu/grch37_tran/Homo_sapiens.GRCh37.75.dna.chromosome.all.fasta --af-of-alleles-not-in-resource 0.000001 --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter --bam-output SC0205_test.recali.m2.bam --output SC0205_test.recali.vcf 
gatk GetPileupSummaries -I Tumor_WES.recali.bam -L /mnt/ngs_data/liu/grch37_tran/Exome_Agilent_V6.bed -V /mnt/ngs_data/liu/grch37_tran/small_exac_common_3_b37.vcf.gz -O pileups.table
gatk SplitReads -I SC0205_test.recali.m2.bam -O /mnt/ngs_data/liu/NDMC/SC0205/ --split-sample
gatk CalculateContamination -I pileups.table -O contamination.table
gatk FilterMutectCalls -R /mnt/ngs_data/liu/grch37_tran/Homo_sapiens.GRCh37.75.dna.chromosome.all.fasta -V SC0205_test.recali.vcf -O filtered.vcf --contamination-table contamination.table
gatk CollectSequencingArtifactMetrics -I SC0205_test.recali.m2.Tumor_WES.bam -O gatk -R /mnt/ngs_data/liu/grch37_tran/Homo_sapiens.GRCh37.75.dna.chromosome.all.fasta --VALIDATION_STRINGENCY=LENIENT
gatk FilterByOrientationBias -V filtered.vcf --output SC0205_test.recali-filtered.vcf -P gatk.pre_adapter_detail_metrics 

