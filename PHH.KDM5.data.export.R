##written to export data for Simon regarding Transport protein genes from fresh samples vs cultured
library(limma)
library(edgeR)
library(Biobase)
library(openxlsx)
####1. PHH, multidonor ####
lg1 <- exprs(readRDS("/fcrbiouatappn01/home/li/Projects/KDM5/rdata/rsem.RDS")$logcpm)
colnames(lg1) <- substring(colnames(lg1),regexpr("/",colnames(lg1))+1)
#meta data
pdat.1 <- pData(readRDS("/fcrbiouatappn01/home/li/Projects/KDM5/rdata/rsem.RDS")$logcpm)
rownames(pdat.1) <- pdat.1$SampleID

d <- matrix(0,nrow=ncol(lg1),ncol=2,dimnames = list(colnames(lg1),c("Fresh","d0")))
d[pdat.1$sample[pdat.1$Day=="Fresh thaw"],"Fresh"] <- 1
d[pdat.1$sample[pdat.1$Day=="d0"&pdat.1$HBV=="none"],"d0"] <- 1
cmat <- matrix(c(-1,1),nrow=2,ncol=1,dimnames = list(colnames(d),c("d0-Fresh")))
fit <- eBayes(contrasts.fit(lmFit(lg1,d),contrasts = cmat))
wb <- createWorkbook("Ricardo")
addWorksheet(wb,"DGE")
dge <- topTable(fit,coef=1,number=Inf)[,c("logFC","P.Value","adj.P.Val")]
writeData(wb,"DGE",dge,rowNames = T)
addWorksheet(wb,"log2CPM")
lcpm <- lg1[,colnames(lg1)%in%pdat.1$sample[pdat.1$HBV=="none"]]
colnames(lcpm) <- substring(colnames(lcpm),1,regexpr("_L",colnames(lcpm))-1)
colnames(lcpm) <- paste(colnames(lcpm),pdat.1[colnames(lcpm),"Day"],pdat.1[colnames(lcpm),"Donor"],sep=".")
writeData(wb,"log2CPM",lcpm,rowNames = T)
saveWorkbook(wb,overwrite = T,
             file="/home/rramirez02/Projects/projects2016.2017/KDM5/fresh.v.plated.PHH.xlsx")

