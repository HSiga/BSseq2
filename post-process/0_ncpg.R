#This code is to calculate number of CPGs with step size from all CPGs bed file of a genome.

library(methylKit)
library(GenomicRanges)
library(GenomicFeatures)
library(data.table)


allCytosines <- read.table("Allcytosines.txt")
allCytosines <-  GRanges(seqnames = allCytosines$V1, ranges = IRanges(start = allCytosines$V2,end=allCytosines$V2), strand = '*') #allCytosines$V3)

num_allCs <- length(allCytosines)
unique(allCytosines@seqnames)
#https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1004868



making_windows_of_nCpGs <- function(single_cpg_list, chunksize = 10, stepsize = 2) {
  
      first_matrix <- unlist(range(equisplit(single_cpg_list, chunksize = chunksize)))

      for (step in (seq(stepsize ,chunksize - stepsize, by = stepsize))) {
        print(step)
        new_step_matrix <- unlist(range(equisplit(tail(single_cpg_list,-step), chunksize = chunksize)))
        first_matrix <- c (first_matrix,new_step_matrix)
        new_step_matrix <- NA
        print(length(first_matrix))
        }
      first_matrix <- sortSeqlevels(first_matrix)
      first_matrix <- sort(first_matrix)
      
      return (first_matrix)
  
}


allCpG_w30_s2 <- making_windows_of_nCpGs(allCytosines, 30, 2)
save(allCpG_w30_s2, file = "allCpG_w30_s2.rdb")



#checks

#test <- making_windows_of_nCpGs(one, 31 ,2)
#head(countOverlaps(test, allCytosines), 100)

#all[countOverlaps(all, allCytosines)==2]
#all[countOverlaps(all, allCytosines)==4]
#all[countOverlaps(all, allCytosines)==6]
#all[countOverlaps(all, allCytosines)==8]
#all[countOverlaps(all, allCytosines)==10]
