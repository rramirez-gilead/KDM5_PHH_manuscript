library(Biobase)
library(limma)
####1. PHH, multidonor EA 15008####
#ex <- readRDS("/home/li/Projects/KDM5/rdata/rsem.RDS")
ex <- readRDS("/fcrbiouatappn01/home/li/Projects/KDM5/rdata/rsem.RDS")
pdat <- pData(ex$logcpm)
ex <- ex$dgelist$counts
colnames(ex) <- substring(colnames(ex),regexpr("/",colnames(ex))+1)
colnames(ex) <- substring(colnames(ex),1,regexpr("_",colnames(ex))-1)
ex <- cbind("Gene"=rownames(ex),ex)
#for(i in 1:ncol(ex)){
#  tmp <- data.frame("Gene"=rownames(ex),
 #                   "Counts"=ex[,i])
  write.table(ex[,c("Gene",as.character(1:53))],file="~/tmp/geo_submission_kdm5/phh.counts.txt",
              quote = F,row.names = F,col.names = T,sep = "\t")
#}
