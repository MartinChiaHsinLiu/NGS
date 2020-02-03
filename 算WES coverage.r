tmp.coverage <- read.table("Tumor_WES.depth.txt",header=F,stringsAsFactors=F)
tmp.coverage <- read.table("Normal_WES.depth.txt",header=F,stringsAsFactors=F)

index.300 <- which(tmp.coverage[,3] < 300)
cov.300 <- tmp.coverage[index.300,3]

index.150 <- which(cov.300 < 150)
cov.150 <- cov.300[index.150]

index.100 <- which(cov.150 < 100)
cov.100 <- cov.150[index.100]

index.50 <- which(cov.100 < 50)
cov.50 <- cov.100[index.50]

index.30 <- which(cov.50 < 30)
cov.30 <- cov.50[index.30]

index.20 <- which(cov.30 < 20)
cov.20 <- cov.30[index.20]

index.10 <- which(cov.20 < 10)
cov.10 <- cov.20[index.10]

cov.0 <- cov.10[which(cov.10 == 0)]

count.300H    <- nrow(tmp.coverage)-length(cov.300)
count.150_300 <- length(cov.300) - length(cov.150)
count.100_150 <- length(cov.150) - length(cov.100)
count.50_100 <- length(cov.100) - length(cov.50)
count.30_50  <- length(cov.50) - length(cov.30)
count.20_30  <- length(cov.30) - length(cov.20)
count.10_20  <- length(cov.20) - length(cov.10)
count.0_10   <- length(cov.10) - length(cov.0)

y<-c()
y[9] <- length(cov.0)/nrow(tmp.coverage)*100
y[8] <- count.0_10   /nrow(tmp.coverage)*100
y[7] <- count.10_20  /nrow(tmp.coverage)*100
y[6] <- count.20_30  /nrow(tmp.coverage)*100 
y[5] <- count.30_50  /nrow(tmp.coverage)*100 
y[4] <- count.50_100 /nrow(tmp.coverage)*100 
y[3] <- count.100_150/nrow(tmp.coverage)*100
y[2] <- count.150_300/nrow(tmp.coverage)*100
y[1] <- count.300H   /nrow(tmp.coverage)*100   
cdf.y <- y[1]
for(i in 2:length(y)){cdf.y[i] <- cdf.y[i-1]+y[i]}

tmp<-data.frame(y)
rownames(tmp)<- rev(c("0","0-10","10-20","20-30","30-50","50-100","100-150","150-300","> 300"))

tmp<-data.frame(cdf.y)
rownames(tmp)<- rev(c("0","> 0","> 10","> 20","> 30","> 50","> 100","> 150","> 300"))








