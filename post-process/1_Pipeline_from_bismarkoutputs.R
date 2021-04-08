#Methylation analysis from BS and oxBS sequencing data.

#Data were processed using TrueMethyl® Data Analysis pipeline provided by CEGX.
#Final outputs from the pre-processed data can be imported to this pipeline in one of the following formats:
#CytosineReports "...CpG_report.txt.gz" or Bismarkcoverage "...bismark.cov.gz" or after conversion
#into Methylkit format  "..methcounts.gz"

#note: Cs in this analysis means CpGs
#put file allCytosines.txt in the same analysis folder (3 columns: chr1 111 +)

#humam.siga@gmail.com

#main packages
library(methylKit)
library(GenomicRanges)
library(GenomicFeatures)
library(data.table)
library(MLML2R)
library(stringr)

setwd(getwd())
input =getwd()
print(input)
#setwd("path.to.files")

#Getting all cpgs in the mouse genome from the Cytosinreport file
#All cpgs in genome
print("Loading all_CpGs file")

#windows of 30 CpGs step 2 by nCpGs.R file
load('allCpG_w30_s2.rdb')
windows_of_cpgs = allCpG_w30_s2

# ---------------------- load bismark outputs for BS and oxBS sequencing files

file_list_BS <- list.files(input, pattern ="^.+_BS_.+gz")
file_list_oxBS <- list.files(input, pattern = "^.+_oxBS_.+gz")

print(data.frame(file_list_BS,file_list_oxBS))

##sample n

for (fileloop in (seq(1,6))) {

  samplename <- word(file_list_BS[fileloop],start = 1, sep =  '_')
  print(paste0(fileloop,"_",samplename))
#}  

single <- c(file_list_BS[fileloop],file_list_oxBS[fileloop])

myobj <- methylKit::methRead(as.list(single), sample.id = list('bs','oxbs'),
                             assembly = "mm10", treatment = c(0,1),
                             context = "CpG", mincov = 0, #will be filtered in unite step
                             pipeline= "bismarkCytosineReport")
#min(myobj[[1]][["coverage"]])
print("Finished loading files")
#Normalize coverage by the samples median coverage
myobj <- methylKit::normalizeCoverage(obj = myobj,method = "median")

#Remove the outlayer coverage
myobj <- filterByCoverage(myobj, lo.count= 0, lo.perc = NULL, hi.count = NULL, hi.perc = 99.9)

#Unite in one DF, (methylBase)
myobj <- methylKit::unite(myobj, destrand = FALSE, min.per.group = 1L, mc.cores = 1) 


#calculate the number of Cs covererd in all samples
# num_cov_Cs <- length(subsetByOverlaps(x = allCytosines, ranges = as(myobj,"GRanges"), minoverlap = 1)) #using myobj


#To re-distribute to single CpG, keep the ones that already had at least one per group
wanted_cpgs <- data.frame(myobj[(myobj$coverage1+myobj$coverage2)>0 ,])[,1:3]
wanted_cpgs <- GRanges(wanted_cpgs)


######################
####### Tiling #######
######################

#method 1 : based on number of CpGs 
tiles <- methylKit::regionCounts(object = myobj,regions = windows_of_cpgs, #introns, using windows of 30 cpg
                                      cov.bases = 0,strand.aware = F)

#Filter low cov in all groups
tiles <- na.omit(tiles)
tiles_100cov <- tiles[tiles$coverage1>=100 & tiles$coverage2>=100 ]

# windows_after <- data.frame(tiles_100cov)[,1:4]
rm(tiles)

########## STEP1 divide the coverage of each window on the remaining CpGS
########## every window should have 30CpGs but some of them have less due to the filtration steps

  #calc total sum
  cc <- countOverlaps( GRanges(tiles_100cov), wanted_cpgs )

  tiles_100cov$coverage1 <- tiles_100cov$coverage1 / cc
  tiles_100cov$numCs1 <-  tiles_100cov$numCs1 / cc
  tiles_100cov$numTs1 <-   tiles_100cov$numTs1 / cc
  tiles_100cov$coverage2 <-  tiles_100cov$coverage2 / cc
  tiles_100cov$numCs2 <-   tiles_100cov$numCs2 / cc
  tiles_100cov$numTs2 <- tiles_100cov$numTs2 / cc
  
  #all overlapped windows merged
  single_cpg <- methylKit::regionCounts(object = tiles_100cov, regions = wanted_cpgs, strand.aware = F)
  

######### STEP2 Re-distibute the numbers in all windows into the CpGs inside it 
######### use only if windows of fix size are used 
######### every CpG is covered with many overlapped windows, thus the mean of them should be calculated

  #calc total sum
  single_cpg <- data.frame(single_cpg)
  #calc number of windows covering every CpG
  tiles_100cov <- as(tiles_100cov,"GRanges")
  wanted_cpg2 <- subsetByOverlaps(wanted_cpgs,ranges = tiles_100cov) #Keep only CpGs that fall in windows.
  seqlevels(wanted_cpg2) <- seqlevelsInUse(wanted_cpg2) #update levels
  
  num_wind <- binnedAverage(bins = wanted_cpg2, numvar = coverage(tiles_100cov), varname = "num_of_wind")$num_of_wind
  num_wind <- num_wind[num_wind!=0]
  #calc average and assign it to meth to estimate 5hmC, you can check length(num_wind) should equal nrow(single_cpg)

  
# single_cpg  
  meth <- cbind(single_cpg[,1:4], single_cpg[,5:10]/num_wind )

# save(meth, file = 'meth')

  rm(single_cpg)


#1 bs
#2 oxbs 
##################################################
# estimation of 5hmc, compaer 1 with 3 

#   eval(parse(text=paste0('meth$numTs',1)))
# 
sampleNo <- samplename

Ts_BS_number <- meth$numTs1
Cs_BS_number <- meth$numCs1
Ts_oxBS_number  <- meth$numTs2
Cs_oxBS_number  <- meth$numCs2

mlml_result <-as.data.frame(
  MLML(U.matrix = Ts_BS_number,  T.matrix = Cs_BS_number,
       L.matrix = Ts_oxBS_number, M.matrix = Cs_oxBS_number))

names(mlml_result) <- paste0(names(mlml_result),'_',sampleNo)

#we need the coverage for every region from the reads
mean_cov <- meth$coverage1 +  meth$coverage2 /2
mean_cov <- data.frame(mean_cov)
names(mean_cov)[1] <- paste0('COVbs_oxbs_sample_',sampleNo) 

assign ( paste0('result_sample_',sampleNo) , cbind.data.frame(meth[,1:4], mean_cov, mlml_result[,1:3]) )

print(head(eval(parse(text= paste0('result_sample_',sampleNo)))))


print("saving file ..")

write.table(x= eval(parse(text= paste0('result_sample_',sampleNo))),
            file = paste0(samplename,"_results_with_cov.txt") , sep = "\t",row.names = F, quote = F)



zipcode <- paste0('gzip ',samplename,"_results_with_cov.txt") 
system(zipcode)

rm(list=setdiff(ls(), c("windows_of_cpgs" ,"file_list_BS" ,"file_list_oxBS")))

}
