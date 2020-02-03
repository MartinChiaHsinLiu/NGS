Sys.setenv(TAR = "/bin/tar")
options(unzip = "internal")

tmp<-read.csv("/home/liu/NDMC/SC0205/output/SNV/MHC_Class_I/Tumor_WES.all_epitopes_flanking.csv",stringsAsFactors=F)
tmp <- tmp[!duplicated(tmp),]


pvacseq_res <- read.delim("/home/liu/NDMC/SC0205/output/SNV/MHC_Class_I/Tumor_WES.all_epitopes.tsv",header=T,stringsAsFactors=F,sep="\t")

write.csv(tmp,"/home/liu/NDMC/SC0205/output/SNV/MHC_Class_I/Tumor_WES.all_epitopes_flanking_parsed.csv",row.names=F,quote=F)

hla<-unique(tmp[,5])
for(i in 1:length(hla)){
    nick.hla <- gsub("\\*","",gsub(":","",gsub("HLA-","",hla[i])))
    tmp.out<-tmp[tmp[,5]==hla[i],]
    write.table(tmp.out,paste0(c("/home/liu/NDMC/SC0205/output/SNV/MHC_Class_I/",nick.hla,".csv"),collapse=""),sep=",",col.names=F,row.names=F,quote=F)
}


python /home/cwlu/EDGE/edgeModelOps.py /home/liu/NDMC/SC0205/output/SNV/MHC_Class_I/A1101.csv /home/liu/NDMC/version1.1/model_A1101_pepdim256_b32_lr3_run0_5_0_.h5
python /home/cwlu/EDGE/edgeModelOps.py /home/liu/NDMC/SC0205/output/SNV/MHC_Class_I/A2402.csv /home/liu/NDMC/version1.1/model_A2402_pepdim256_b32_lr3_run0_5_0_.h5
python /home/cwlu/EDGE/edgeModelOps.py /home/liu/NDMC/SC0205/output/SNV/MHC_Class_I/B1501.csv /home/liu/NDMC/version1.1/model_B1501_pepdim256_b32_lr3_run0_2_0_.h5
python /home/cwlu/EDGE/edgeModelOps.py /home/liu/NDMC/SC0205/output/SNV/MHC_Class_I/B5801.csv /home/liu/NDMC/version1.1/model_B5801_pepdim256_b32_lr3_run0_10_0_.h5
python /home/cwlu/EDGE/edgeModelOps.py /home/liu/NDMC/SC0205/output/rescue/MHC_Class_I/A1101.csv /home/liu/NDMC/version1.1/model_A1101_pepdim256_b32_lr3_run0_5_0_.h5
python /home/cwlu/EDGE/edgeModelOps.py /home/liu/NDMC/SC0205/output/rescue/MHC_Class_I/A2402.csv /home/liu/NDMC/version1.1/model_A2402_pepdim256_b32_lr3_run0_5_0_.h5
python /home/cwlu/EDGE/edgeModelOps.py /home/liu/NDMC/SC0205/output/rescue/MHC_Class_I/B1501.csv /home/liu/NDMC/version1.1/model_B1501_pepdim256_b32_lr3_run0_2_0_.h5
python /home/cwlu/EDGE/edgeModelOps.py /home/liu/NDMC/SC0205/output/rescue/MHC_Class_I/B5801.csv /home/liu/NDMC/version1.1/model_B5801_pepdim256_b32_lr3_run0_10_0_.h5

tmp<-read.csv("/home/liu/NDMC/SC0205/output/SNV/MHC_Class_I/Tumor_WES.all_epitopes_flanking.csv",stringsAsFactors=F)
tmp <- tmp[!duplicated(tmp),]
#tmp[tmp[,3]=="X",3]<-"----------"

pvacseq_res <- read.delim("/home/liu/NDMC/SC0205/output/SNV/MHC_Class_I/Tumor_WES.all_epitopes.tsv",header=T,stringsAsFactors=F,sep="\t")

res.files <- list.files(pattern="predictResult.csv")
for(i in 1:length(res.files)){
res.tmp <- read.csv(res.files[i],header=F,stringsAsFactors=F)
if(i==1){res<-res.tmp}else{res<-rbind(res,res.tmp)}
}
res<-res[,2:8]
res[,2]<-gsub("Z","-",res[,2])

colnames(res) <- c("MT.Epitope.Seq","flanking","HLA.Allele","protein","mono","log10TypescriptExpression","score")

res.score <- merge(tmp,res)

write.csv(res.score,"/home/liu/NDMC/SC0205/output/SNV/MHC_Class_I/Tumor_WES.all_epitopes_flanking_presentationscore.csv",row.names=F,quote=F)

res2 <- res[,c(1,2,3,6,7)]
colnames(res2)[4]<-"Transcript.Expression"
res2[,4] <- 10^res2[,4]
pvacseq_res.score <- merge(pvacseq_res,res2)

write.table(pvacseq_res.score ,"/home/liu/NDMC/SC0205/output/SNV/MHC_Class_I/Tumor_WES.all_epitopes_score.tsv",sep="\t",row.names=F,quote=F)









export AG_DATA_DIR="/home/liu/antigen.garnish"
export PATH=$PATH:/home/liu/antigen.garnish/ncbi-blast-2.7.1+/bin
tmp <- read.csv("/home/liu/EDGE/CJWu_2017.csv",header=T,stringsAsFactors=F)
res.diss <- garnish_dissimilarity(v,db = "human")
res.iedb <- iedb_score(v,db = "human")
res.diss@.Data[[1]]
res.diss@.Data[[2]]


library(magrittr)
library(data.table)
library(antigen.garnish)
source("/home/liu/garnish_dissimilarity.r") #因相容性問題改過data.table::xxxxx => xxxxx 
v <- read.csv("/home/liu/EDGE/immunognenicity/gritstone_response_input.csv",header=F,stringsAsFactors=F)[,1]
res.diss <- garnish_dissimilarity(v,db = "human")
res.iedb <- iedb_score(v,db = "human")


‘Rsamtools’,  ‘VGAM’, ‘GenomicAlignments’








