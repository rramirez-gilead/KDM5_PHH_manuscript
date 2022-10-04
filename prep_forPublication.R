library(Biobase)
library(limma)
library(GSVA)
####1. PHH, multidonor EA 15008####
#import data from Li
ex <- readRDS("/fcrbiouatappn01/home/li/Projects/KDM5/rdata/rsem.RDS")
pdat <- pData(ex$logcpm)
ex <- ex$dgelist$counts
colnames(ex) <- substring(colnames(ex),regexpr("/",colnames(ex))+1)
colnames(ex) <- substring(colnames(ex),1,regexpr("_",colnames(ex))-1)
ex <- cbind("Gene"=rownames(ex),ex)
#run once
if(interactive()){
  write.table(ex[,c("Gene",as.character(1:53))],file="~/tmp/geo_submission_kdm5/phh.counts.txt",
              quote = F,row.names = F,col.names = T,sep = "\t")
  #load fpkm data
  fpkm <- read.table("/fcrbiouatappn01/home/li/Projects/KDM5/Batch.2015-02-27/rsem_gfpkm_matrix.txt")
  loadGeneSets <- function(infile="h.all.v7.4.symbols.gmt",
                           indir="/fcrbiouatappn01/home/rramirez02/Projects/HBsInhibitor/"){
    sets <- alist()
    sfile <-  file(paste(indir,infile,sep="/"),open = "r")
    while(length(l <- readLines(sfile,n=1)) > 0){
      tmp <- strsplit(l,split="\t")[[1]]
      sets[[tmp[1]]] <- tmp[3:length(tmp)]
    }
    close(sfile)
    rm(tmp)
    return(sets)
  }
  
  genesets <- list("hallmark"=loadGeneSets(),
                   "c2"=loadGeneSets("c2.all.v7.4.symbols.gmt"),
                   "c3"=loadGeneSets("c3.all.v7.4.symbols.gmt"),
                   "c5"=loadGeneSets("c5.all.v7.4.symbols.gmt"),
                   "tft"=loadGeneSets("c3.tft.v7.4.symbols.gmt"))
  
  rungsva <- function(set="hallmark",indata=fpkm,allgenes=unique(unlist(genesets))){
    res <- gsva(as.matrix(fpkm),genesets[[set]])
    return(res)
  }
  gsv <- lapply(names(genesets),rungsva);names(gsv) <- names(genesets)
  saveRDS(gsv,"gsva.RDS")
}else{
  gsv <- readRDS("gsva.RDS")
}

##gene set enrichment analysis


runlim <- function(indata=lcpm,minmx=NULL,trnd=F){
  # if(!is.null(minmx)){#filter genes with low expression
  #   #calculate cohort means
  #   means <- matrix(0,nrow=nrow(indata),ncol=4,
  #                   dimnames = list(rownames(lcpm),unique(cohorts)))
  #   for(r in rownames(means)){
  #     means[r,] <- sapply(colnames(means),function(x) mean(indata[r,names(cohorts)[cohorts%in%x]]))
  #   }
  #   keep <- apply(means,1,max)>minmx
  #   indata <- indata[keep,]
  # }
  # d <- matrix(0,nrow=ncol(indata),ncol=4,
  #             dimnames = list(colnames(indata),
  #                             unique(cohorts)))
  # for(s in rownames(d)){
  #   d[s,cohorts[s]] <- 1
  # }
  # cmat <- matrix(0,nrow=ncol(d),ncol=3,
  #                dimnames = list(colnames(d),c("GS353-veh","GS635-ABT","ABT-veh")))
  # cmat[c(2,8,11)] <- 1
  # cmat[c(1,7,9)] <- -1
  # fit <- lmFit(indata,d)
  # res <- eBayes(contrasts.fit(fit,contrasts = cmat),trend = trnd)
  # attributes(res)$fit <- fit
  # return(res)
}
lim <- list("genes"=runlim(minmx=1,trnd=T),
            "hallmark"=runlim(gsv$hallmark),
            "c2"=runlim(gsv$c2),"c3"=runlim(gsv$c3),
            "c5"=runlim(gsv$c5),
            "tft"=runlim(gsv$tft),
            "kegg"=runlim(gsv$c2[grep("^KEGG",rownames(gsv$c2)),]),
            "wp"=runlim(gsv$c2[grep("^WP",rownames(gsv$c2)),]))