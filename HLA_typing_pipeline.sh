#Install Opitype 與必要套件:
#optitype用bioconda安裝
#conda install -c bioconda/label/cf201901 -c bioconda -c conda-forge optitype
#安裝完以後要把/mnt/ngs_data/liu/.conda/pkgs/optitype-1.2-py35_0/share/optitype-1.2-0/config.ini裡面的
#razers3=/path/to/razers3改成razers3=razers3
#Argument說明:
#-i 是input files (這裡要用的是normal part的部份)
#-o 是output dir
conda activate optitype
cd /mnt/ngs_data/liu/NDMC/SC0205/

python /home/ubuntu/anaconda3/pkgs/optitype-1.3.2-py27_1/bin/OptiTypePipeline.py -i Normal_WES_R1.paired.trim.fastq.gz Normal_WES_R2.paired.trim.fastq.gz --dna -v -o /mnt/ngs_data/liu/NDMC/SC0205/opitype/



