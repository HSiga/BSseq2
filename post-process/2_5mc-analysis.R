# Calculate DMRs for 5mC

#module load bioinfo-tools
#module load R/3.4.3
#module load R_packages/3.4.3
library("bsseq")
setwd("Path/to/cov.gz") # from previous step
 
                                                   #5MC analysis
                                                   #READ DATA
 #trace(bsseq:::read.bismarkCovRaw, edit = TRUE)
 
 #Adjust the function for importing the files
 read.bismarkCovRaw <- function (thisfile, thisSampleName, rmZeroCov, BACKEND = NULL) 
 {
   if (isGzipped(thisfile)) {
     thisfile <- gunzip(thisfile, temporary = TRUE, overwrite = TRUE, remove = FALSE)
   }
   out <- fread(thisfile)
   out <- na.omit(out)
   if (ncol(out) != 8L) {
     stop("unknown file format")
   }
   gr <- GRanges(seqnames = out[[1L]], ranges = IRanges(start = out[[2L]], width = 1L))
   M <- as.matrix(out[[6L]] * out[[5L]])
   Cov <- as.matrix(out[[5L]] )
   M <- realize(M, BACKEND = BACKEND)
   Cov <- realize(Cov, BACKEND = BACKEND)
   BSseq(gr = gr, M = M, Cov = Cov, sampleNames = thisSampleName, 
         rmZeroCov = rmZeroCov)
 }
 
 #To reassign the loading function
 tmpfun <- get('read.bismarkCovRaw', envir = asNamespace('bsseq'))
 environment(read.bismarkCovRaw) <- environment(tmpfun)
 assignInNamespace('read.bismarkCovRaw', read.bismarkCovRaw, ns = 'bsseq')
 
 #load files
 bismarkBSseq <- read.bismark(files=c( 'WT1_results_with_cov.txt.gz',
                                       'WT2_results_with_cov.txt.gz',
                                       'WT3_results_with_cov.txt.gz',
                                       'HOM1_results_with_cov.txt.gz',
                                       'HOM2_results_with_cov.txt.gz',
                                       'HOM3_results_with_cov.txt.gz'),
                              sampleNames =c("WT2","WT3","WT3","HOM1","HOM2","HOM3"),
                              rmZeroCov=TRUE,strandCollapse=FALSE,fileType="cov",mc.cores=15,verbose=TRUE)
 
 #DEFINE GROUPS
 pData(bismarkBSseq)$Type<-c("WT","WT","WT","HOM","HOM","HOM")
 bismarkBSseq
 
 #An object of type 'BSseq' with
 #  40357584 methylation loci
 #  6 samples
 #has not been smoothed
 #All assays are in-memory
 
 #CALCULATE MEAN COVERAGE PER SAMPLE
 print(round(colMeans(getCoverage(bismarkBSseq)),1))
 # WT2     WT3    WT3B HOM2355 HOM2356 HOM2357 
 # 14.8    16.9    15.5    14.3    13.8    14.3 
 
 
 print(paste0("TOTAL NUMBER OF CPGs: ",length(bismarkBSseq)))
 #"TOTAL NUMBER OF CPGs: 40357584"
 print(paste0("NUMBER OF CPGs COVERED BY AT LEAST TWO READS IN ALL SAMPLES: ",
              sum(rowSums(getCoverage(bismarkBSseq)>=2)==6)))
 #[1] "NUMBER OF CPGs COVERED BY AT LEAST TWO READS IN ALL SAMPLES: 34 947 121"
 print(paste0("NUMBER OF CPGs WITH ZERO COVERAGE IN ALL SAMPLES: ",sum(rowSums(getCoverage(bismarkBSseq))==0)))

  #load('5mC_bismarkBSseq.rdb')
  pdf("5mC_output_dmrs_3WT_vs_3HOM.pdf",paper="a4r",height=210,width=297)
  dmrs<-list()
  dmrs0<-list()
  bismarkBSseq.fit<-list()
  bismarkBSseq.fit.cov<-list()
  bismarkBSseq.fit.tstat<-list()
  
#save (bismarkBSseq, file = '5mC_bismarkBSseq.rdb')

for(i in c(1:19,'X'))
{
  j <- paste0('chr',i)
  gc()
  print(paste0("WORKING WITH CHROMOSOME ",i))
  bismarkBSseq.subset<-chrSelectBSseq(bismarkBSseq,seqnames=as.character(j),order=TRUE)
  if (length(bismarkBSseq.subset) == 0) next #skip chrs where there is no hits
  bismarkBSseq.fit[[i]]<-BSmooth(bismarkBSseq.subset,mc.cores=15,verbose=TRUE)
  
  print("CALCULATE COVERAGE AND SELECT ONLY WELL COVERED SITES")
  bismarkBSseq.fit.cov[[i]]<-getCoverage(bismarkBSseq.fit[[i]])
  
  #KEEP ONLY LOCI WHERE BOTH INFANT SAMPLES AND BOTH NON-INFANT SAMPLES HAVE AT LEAST 2X COVERAGE
  keepLoci<-which(rowSums(bismarkBSseq.fit.cov[[i]][,bismarkBSseq.fit[[i]]$Type=="WT"]>=2)>=3 &
                    rowSums(bismarkBSseq.fit.cov[[i]][,bismarkBSseq.fit[[i]]$Type=="HOM"]>=2)>=3)
  bismarkBSseq.fit[[i]]<-bismarkBSseq.fit[[i]][keepLoci,]

  if (length(bismarkBSseq.fit[[i]]) == 0) next #skip chrs where there is no hits

  print("CALCULATE T-STAT")
  bismarkBSseq.fit.tstat[[i]] <- BSmooth.tstat(bismarkBSseq.fit[[i]],
                                               group1 = c("WT2","WT3","WT3B"),
                                               group2 = c("HOM2355","HOM2356","HOM2357"),
                                               verbose = TRUE)
  plot(bismarkBSseq.fit.tstat[[i]],main=paste0("CHROMOSOME ",i))
  
  print("CALCULATE DMRs")
  dmrs0[[i]]<-dmrFinder(bismarkBSseq.fit.tstat[[i]])
  
  #FILTER OUT DMRs THAT HAVE LESS THAN 3 CPGs AND WT VS. HOM DIFFERENCE BELOW 10%
  dmrs[[i]]<-subset(dmrs0[[i]], n>=3 & abs(meanDiff)>=0.1)
  print(paste0("FOUND ",nrow(dmrs[[i]])," DMRs FOR CHROMOSOME ",i))
  print(head(dmrs[[i]]))
  
  pData<-pData(bismarkBSseq.fit[[i]])
  pData$col<-c("blue","blue","blue","red","red","red")
  pData(bismarkBSseq.fit[[i]])<-pData
  plotRegion(bismarkBSseq.fit[[i]],dmrs[[i]][1,],extend=5000,addRegions=dmrs[[i]],main=paste0("CHROMOSOME ",i))
  legend("left",inset=.02,col=c("blue","red"),lty=c(1,1),c("WT","HOM"))
  
  pdf(file=paste0("5mC_dmrs_top100_3WT_vs_3HOM_",i,".pdf"),paper="a4r",height=210,width=297)
  plotManyRegions(bismarkBSseq.fit[[i]],dmrs[[i]][1:100,],extend=5000,
                  addRegions=dmrs[[i]],main=paste0("CHROMOSOME ",i))
  legend("left",inset=.02,col=c("blue","red"),lty=c(1,1),c("WT","HOM"))
  dev.off()
  
  }
dev.off()


output<-Reduce(rbind,dmrs)

#RANK DMRs WITH SUM OF T-STATISTICS AT EACH CPG WEIGHTED BY NUMBER OF CPGs
output<-output[order(-abs(output$areaStat)),]
print(head(output,10))
write.table(output,file="5mC_output_dmrs_3WT_vs_3HOM.txt",
            col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")


