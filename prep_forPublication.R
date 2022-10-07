library(Biobase)
library(limma)
library(GSVA)
#make figure or of hbv RNA from RNA-Seq
##send to Uli as TIFF
####1. PHH, multidonor EA 15008####
#import data from Li
ex <- readRDS("/fcrbiouatappn01/home/li/Projects/KDM5/rdata/rsem.RDS")
pdat <- pData(ex$logcpm)
ex <- ex$dgelist$counts
colnames(ex) <- substring(colnames(ex),regexpr("/",colnames(ex))+1)
colnames(ex) <- substring(colnames(ex),1,regexpr("_",colnames(ex))-1)
ex <- cbind("Gene"=rownames(ex),ex)
#run once
if(!interactive()){
  write.table(ex[,c("Gene",as.character(1:53))],file="/fcrbiouatappn01/home/rramirez02/tmp/geo_submission_kdm5/phh.counts.txt",
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
  
  ##align to HBV####
  runHBValign <- function(id="7",
                          fqdir="/fcrbiouatappn01/home/li/Projects/KDM5/Samples.2015-02-27/20150310-H42010.FASTQs"){
    fq1 <- grep("1.clipped.fastq.gz",dir(fqdir,pattern=paste(id,"_",sep="")),value=T)
    fq1 <- paste(fqdir,fq1,sep="/")
    fq2 <- sub("1.clipped","2.clipped",fq1,fixed=T)
    outsam <- paste("sam/",id,".hbv.sam",sep="")
    cmd <- paste("/export/home/rramirez02/home_rramirez02/bin/bwa mem -t 6",
                 "/export/home/rramirez02/home_rramirez02/Virus/HBV.AD38.fa",
                 fq1,fq2,"| grep HBV >",outsam)
    system(cmd)
    #outbam <- sub("sam","bam",outsam)
    #cmd <- paste("/export/home/rramirez02/home_rramirez02/bin/samtools sort",)
  }
  lidir <- "/fcrbiouatappn01/home/li/Projects/KDM5/Samples.2015-02-27/"
  fqdirs <- dir(lidir,pattern="FASTQs")
  for(i in fqdirs){
    ids <- dir(paste(lidir,i,sep=""))
    ids <- sort(unique(substring(ids,1,regexpr("_",ids)-1)))
    for(tmp in ids){
      runHBValign(tmp,paste(lidir,i,sep=""))
    }
  }
  #create HBV counts
  samdir <- "sam/"
  sams <- dir(samdir,pattern="sam")
  hbvcounts <- NULL#will store raw counts
  for(s in sams){
    cmd <- paste("cut -f 1 ",samdir,s," | grep -v @ | uniq | sort | uniq | wc -l", sep="")
    hbvcounts[s] <- as.numeric(system(cmd,intern = T))
  }
  hbvcounts <- as.numeric(hbvcounts)
  names(hbvcounts) <- sub(".hbv.sam","",sams,fixed = T)
  #get library sizes
  raw <- read.table("/fcrbiouatappn01/home/li/Projects/KDM5/Batch.2015-02-27/rsem_gcounts_matrix.txt",
                    sep="\t",header=T,row.names = 1)
  libsizes <- colSums(raw)
  hbvlog2 <- NULL
  for(i in names(libsizes)){
    tmp <- substring(i,regexpr(".",i,fixed=T)+1,regexpr("_",i,fixed=T)-1)
    hbvlog2[tmp] <- log2((.25+hbvcounts[tmp])/(libsizes[i]/1000000))
  }
  hbvcounts <- list("rawcounts"=hbvcounts,
                    "libsizes"=libsizes,
                    "log2CPM"=hbvlog2)
  saveRDS(hbvcounts,"hbv.RDS")
}else{
  gsv <- readRDS("gsva.RDS")
  hbvcounts <- readRDS("hbv.RDS")
}

##gene set enrichment analysis####

runlim <- function(indata=gsv$hallmark,trnd=F){
  res <- alist()
  for(day in c("d1","d3","d10","d13")){
    dpd <- subset(pdat,Day==day)
    d <- matrix(0,nrow=nrow(dpd),ncol=6,
                dimnames = list(rownames(dpd),
                                c("d8181","d4239",
                                  paste("dose",c(0,0.03,0.3,10),sep = ""))))
    for(don in c("8181","4239")){
      d[rownames(dpd)[dpd$Donor==don],paste("d",don,sep="")] <- 1
    }
    for(dose in c(0,0.03,0.3,10)){
      d[rownames(dpd)[dpd$Dose==dose],paste("dose",dose,sep="")] <- 1
    }
    rownames(d) <- sub("/",".",rownames(d))
    cmat <- matrix(0,nrow=6,ncol=3,
                   dimnames=list(colnames(d),
                                 c("0.03-ctrl","0.3-ctrl","10-ctrl")))
    cmat["dose0",] <- -1
    for(dose in c(0.3,0.03,10)){
      cmat[paste("dose",dose,sep=""),paste(dose,"-ctrl",sep="")] <- 1
    }
   fit <- lmFit(indata[,rownames(d)],d)
   cfit <- contrasts.fit(fit,cmat)
   efit <- eBayes(cfit,trend=trnd)
   attributes(efit)$fit <- fit
   attributes(efit)$cfit <- cfit
   res[[day]] <- efit
  }
  return(res)
}
lim <- list(#"genes"=runlim(minmx=1,trnd=T),
            "hallmark"=runlim(gsv$hallmark),
            "c2"=runlim(gsv$c2),"c3"=runlim(gsv$c3),
            "c5"=runlim(gsv$c5),
            "tft"=runlim(gsv$tft),
            "kegg"=runlim(gsv$c2[grep("^KEGG",rownames(gsv$c2)),]),
            "wp"=runlim(gsv$c2[grep("^WP",rownames(gsv$c2)),]))


pheatgsv <- function(inlim=lim$hallmark,mask="HALLMARK_|KEGG_|WP_",showp=T,adjustP=T,filterF=F,...){
  require(pheatmap)
  coefs <- NULL
  for(day in names(inlim)){
    tmp <- inlim[[day]]$coefficients
    colnames(tmp) <- paste(day,colnames(tmp))
    coefs <- cbind(coefs,tmp)
  }
  acol <- data.frame( "HBV"=NA,
    "Day"=substring(colnames(coefs),1,regexpr("[ ]",colnames(coefs))-1),
                     "Dose"=sub("-ctrl","",substring(colnames(coefs),regexpr("[ ]",colnames(coefs))+1)),
                     row.names = colnames(coefs))
  acol$Day <- factor(acol$Day,levels=c("d1","d3","d10","d13"))
  acol <- acol[order(acol$Day),]
  for(i in 1:nrow(acol)){
    ctrlid <- pdat$SampleID[as.character(pdat$Day)==as.character(acol[i,"Day"]) &
                              pdat$Dose==0]
    treatid <- pdat$SampleID[as.character(pdat$Day)==as.character(acol[i,"Day"]) &
                        as.character(pdat$Dose)==acol[i,"Dose"]]
    #for(d in doses){
      #treat <- pdat$SampleID[as.character(pdat$Day)==as.character(acol[i,"Day"]) &
      #                         pdat$Dose==d]
      acol[i,"HBV"] <- mean(hbvcounts$log2CPM[as.character(treatid)])-mean(hbvcounts$log2CPM[as.character(ctrlid)])
    #}
  }
  acol <- acol[,c("HBV","Dose")]
  colnames(acol) <- c("logFC HBV","Dose")
  mx <- max(abs(coefs))
  brks <- (mx*(-32:32))/32
  if(showp){
    ps<- NULL
    for(day in names(inlim)){
      tmp <- inlim[[day]]$p.value
      colnames(tmp) <- paste(day,colnames(tmp))
      if(adjustP){
        for(dose in colnames(tmp)){
          tmp[,dose] <- p.adjust(tmp[,dose])
        }
      }
      ps <- cbind(ps,tmp)
    }
    sig <- ps
    sig[ps<=0.05] <- "*"
    sig[ps>0.05] <- ""
    ps <- sig
  }else{
    ps <- F

  }
  
  rownames(coefs) <-sub(mask,"",rownames(coefs))
  rownames(coefs) <- gsub("_"," ",rownames(coefs))

  pheatmap(coefs[,rownames(acol)],
           color = colorRampPalette(c("navy","blue","gray95","red","darkred"))(64),
           breaks = brks,
           display_numbers = ps,
          border_color = NA,
          gaps_col = c(3,6,9),
          annotation_colors = list("Day"=c("d1"="#10c070","d3"="#058045","d10"="#056035","d13"="#054020"),
                                   "Dose"=c("0.03"="gray80","0.3"="gray40","10"="gray20"),
                                   "logFC HBV"=c("#00d000","gray10")),
           annotation_col = acol,cluster_cols = F,show_colnames = F,...)
}
pheatgsv(fontsize_row=7,fontsize=7,fontsize_number=10,filename="hallmark.tiff")
#pheatgsv(lim$wp,fontsize_row=6,show_rownames=F,fontsize=7,showp = F)
#pheatgsv(lim$kegg,fontsize_row=6,show_rownames=F,fontsize=7,showp = F)
#pheatgsv(lim$kegg,fontsize_row=6,show_rownames=F,fontsize=7,showp = T)
