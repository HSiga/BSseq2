setwd("/proj/snic2021-22-81/nobackup/mC_methyl_regions_CpG")
diff_meth_sites<-readLines("/proj/snic2021-22-81/nobackup/mc_CpG_sites_10perc_diff.txt")
func_elements<-c("CDS","exon","intron","promoter","enhancer","CTCF","TF","open_chrom")
for(i in 1:length(func_elements))
{
  df<-read.delim(paste0("df_5mc_CpG_methyl_",func_elements[i],".txt"),header=TRUE,row.names=1,check.names=FALSE,sep="\t")
  df<-df[rownames(df)%in%diff_meth_sites,]
  write.table(df,file=paste0("/proj/snic2021-22-81/nobackup/mC_methyl_regions_diff_meth_CpG/df_5mc_CpG_methyl_",func_elements[i],".txt"),
              col.names=TRUE,row.names=TRUE,quote=FALSE,sep="\t")
  print(paste0("Finished ",func_elements[i]))
}



setwd("/proj/snic2021-22-81/nobackup/hmC_methyl_regions_CpG")
diff_meth_sites<-readLines("/proj/snic2021-22-81/nobackup/hmc_CpG_sites_10perc_diff.txt")
func_elements<-c("CDS","exon","intron","promoter","enhancer","CTCF","TF","open_chrom")
for(i in 1:length(func_elements))
{
  df<-read.delim(paste0("df_5hmc_CpG_methyl_",func_elements[i],".txt"),header=TRUE,row.names=1,check.names=FALSE,sep="\t")
  df<-df[rownames(df)%in%diff_meth_sites,]
  write.table(df,file=paste0("/proj/snic2021-22-81/nobackup/hmC_methyl_regions_diff_meth_CpG/df_5hmc_CpG_methyl_",func_elements[i],".txt"),
              col.names=TRUE,row.names=TRUE,quote=FALSE,sep="\t")
  print(paste0("Finished ",func_elements[i]))
}

