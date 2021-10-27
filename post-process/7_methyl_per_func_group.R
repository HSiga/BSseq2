############################ BUILDING METHYLATION LEVELS BASED ON DMR WITHIN FUNCTIONAL REGIONS ########################################
mc_matrix<-read.delim("/home/nikolay/WABI/Misc/Zaher/mC_matrix.txt",header=TRUE,sep="\t")
mc_coords<-matrix(unlist(strsplit(rownames(mc_matrix),"_")),ncol=3,byrow=TRUE)
mc_matrix$chr<-mc_coords[,1]
mc_matrix$start<-mc_coords[,2]
mc_matrix$end<-mc_coords[,3]

annot<-c("CDS","exon","intron","start_codon","stop_codon","promoter","enhancer","CTCF","TF","open_chrom")
#annot<-c("start_codon","stop_codon","promoter","enhancer","CTCF","TF","open_chrom")
#annot<-c("TF","open_chrom")
for(k in annot)
{
print(paste0("WORKING WITH: ",k," ######################"))
regions<-read.delim(paste0("/home/nikolay/WABI/Misc/Zaher/reference_data/df_",k,".txt"),header=FALSE,sep="\t")
random_indices<-sample(1:dim(regions)[1],1000)
df_region_meth<-data.frame(matrix(ncol=dim(mc_matrix)[2]))
colnames(df_region_meth)<-colnames(mc_matrix)
df_region_meth<-df_region_meth[-1,]
j<-1
for(i in random_indices)
{
region_meth<-mc_matrix[as.character(mc_matrix$chr)==as.character(regions[i,]$V1) & as.numeric(as.character(mc_matrix$start))>=as.numeric(as.character(regions[i,]$V4)) & as.numeric(as.character(mc_matrix$end))<=as.numeric(as.character(regions[i,]$V5)),]
region_meth$chr<-NULL
region_meth$start<-NULL
region_meth$end<-NULL
region_mean_meth<-colMeans(region_meth)
region_mean_meth<-as.data.frame(region_mean_meth)
colnames(region_mean_meth)<-paste0(as.character(regions[i,]$V1),"_",as.character(regions[i,]$V4),"_",as.character(regions[i,]$V5))
region_mean_meth<-t(region_mean_meth)
df_region_meth<-rbind(df_region_meth,region_mean_meth)
print(paste0("FINISHED ",j," REGIONS"))
j<-j+1
}
write.table(df_region_meth,file=paste0("/home/nikolay/WABI/Misc/Zaher/mC_methyl_regions/df_5mc_methyl_",k,".txt"),col.names=TRUE,row.names=TRUE,quote=FALSE,sep="\t")
}


#################### PLOT HOM VS WT METHYLATION LEVELS ###############################################
mean_HOM_mc_reg<-vector()
mean_WT_mc_reg<-vector()
sd_HOM_mc_reg<-vector()
sd_WT_mc_reg<-vector()
annot<-c("CDS","exon","intron","promoter","enhancer","CTCF","TF","open_chrom")
for(k in annot)
{
df<-read.delim(paste0("/home/nikolay/WABI/Misc/Zaher/mC_methyl_regions/df_5mc_methyl_",k,".txt"),header=TRUE,sep="\t")
df_HOM<-subset(df,select=colnames(df)[grepl("HOM",colnames(df))])
df_WT<-subset(df,select=colnames(df)[grepl("WT",colnames(df))])
mean_HOM_mc_reg<-append(mean_HOM_mc_reg,mean(as.matrix(df_HOM),na.rm=TRUE))
mean_WT_mc_reg<-append(mean_WT_mc_reg,mean(as.matrix(df_WT),na.rm=TRUE))
sd_HOM_mc_reg<-append(sd_HOM_mc_reg,sd(as.matrix(df_HOM),na.rm=TRUE))
sd_WT_mc_reg<-append(sd_WT_mc_reg,sd(as.matrix(df_WT),na.rm=TRUE))
print(paste0("FINISHED ",k))
}

mean_HOM_hmc_reg<-vector()
mean_WT_hmc_reg<-vector()
sd_HOM_hmc_reg<-vector()
sd_WT_hmc_reg<-vector()
annot<-c("CDS","exon","intron","promoter","enhancer","CTCF","TF","open_chrom")
for(k in annot)
{
df<-read.delim(paste0("/home/nikolay/WABI/Misc/Zaher/hmC_methyl_regions/df_5hmc_methyl_",k,".txt"),header=TRUE,sep="\t")
df_HOM<-subset(df,select=colnames(df)[grepl("HOM",colnames(df))])
df_WT<-subset(df,select=colnames(df)[grepl("WT",colnames(df))])
mean_HOM_hmc_reg<-append(mean_HOM_hmc_reg,mean(as.matrix(df_HOM),na.rm=TRUE))
mean_WT_hmc_reg<-append(mean_WT_hmc_reg,mean(as.matrix(df_WT),na.rm=TRUE))
sd_HOM_hmc_reg<-append(sd_HOM_hmc_reg,sd(as.matrix(df_HOM),na.rm=TRUE))
sd_WT_hmc_reg<-append(sd_WT_hmc_reg,sd(as.matrix(df_WT),na.rm=TRUE))
print(paste0("FINISHED ",k))
}

my_barplot<-rbind(mean_HOM_mc_reg,mean_WT_mc_reg)
barplot(my_barplot,beside=T,col=c("blue","orange"),ylab="Methylation Fraction",names=annot,ylim=c(0,1),main="5mC")
legend("topright",c("HOM","WT"),inset=0.02,fill=c("blue","orange"),cex=1.2)


my_barplot<-rbind(mean_HOM_hmc_reg,mean_WT_hmc_reg)
barplot(my_barplot,beside=T,col=c("blue","orange"),ylab="Methylation Fraction",names=annot,ylim=c(0,0.1),main="5hmC")
legend("topright",c("HOM","WT"),inset=0.02,fill=c("blue","orange"),cex=1.2)

#library("stringr")
#mc_matrix<-mc_matrix[str_count(rownames(mc_matrix),"_")==2,]
#mc_matrix<-na.omit(mc_matrix)

#cut -f1 mC_matrix.txt | sed '1d' > mc_coords.txt
#mc_coords<-readLines("mc_coords.txt")
#mc_coords<-mc_coords[str_count(mc_coords,"_")==2]
#mc_coords<-matrix(unlist(strsplit(mc_coords,"_")),ncol=3,byrow=TRUE)









exon_coords<-read.delim("/proj/snic2021-22-81/nobackup/CpG_func_reg_annotated/df_exon_coords.txt",header=FALSE,sep="\t")
mc_exon<-mc[rownames(mc)%in%paste0(exon_coords$V1,"_",exon_coords$V2,"_",exon_coords$V3),]
hmc_exon<-hmc[rownames(hmc)%in%paste0(exon_coords$V1,"_",exon_coords$V2,"_",exon_coords$V3),]
write.table(mc_exon,file="df_5mc_CpG_methyl_exon.txt",col.names=TRUE,row.names=TRUE,quote=FALSE,sep="\t")
write.table(hmc_exon,file="df_5hmc_CpG_methyl_exon.txt",col.names=TRUE,row.names=TRUE,quote=FALSE,sep="\t")

promoter_coords<-read.delim("/proj/snic2021-22-81/nobackup/CpG_func_reg_annotated/df_promoter_coords.txt",header=FALSE,sep="\t")
mc_promoter<-mc[rownames(mc)%in%paste0(promoter_coords$V1,"_",promoter_coords$V2,"_",promoter_coords$V3),]
hmc_promoter<-hmc[rownames(hmc)%in%paste0(promoter_coords$V1,"_",promoter_coords$V2,"_",promoter_coords$V3),]
write.table(mc_promoter,file="df_5mc_CpG_methyl_promoter.txt",col.names=TRUE,row.names=TRUE,quote=FALSE,sep="\t")
write.table(hmc_promoter,file="df_5hmc_CpG_methyl_promoter.txt",col.names=TRUE,row.names=TRUE,quote=FALSE,sep="\t")

enhancer_coords<-read.delim("/proj/snic2021-22-81/nobackup/CpG_func_reg_annotated/df_enhancer_coords.txt",header=FALSE,sep="\t")
mc_enhancer<-mc[rownames(mc)%in%paste0(enhancer_coords$V1,"_",enhancer_coords$V2,"_",enhancer_coords$V3),]
hmc_enhancer<-hmc[rownames(hmc)%in%paste0(enhancer_coords$V1,"_",enhancer_coords$V2,"_",enhancer_coords$V3),]
write.table(mc_enhancer,file="df_5mc_CpG_methyl_enhancer.txt",col.names=TRUE,row.names=TRUE,quote=FALSE,sep="\t")
write.table(hmc_enhancer,file="df_5hmc_CpG_methyl_enhancer.txt",col.names=TRUE,row.names=TRUE,quote=FALSE,sep="\t")

CTCF_coords<-read.delim("/proj/snic2021-22-81/nobackup/CpG_func_reg_annotated/df_CTCF_coords.txt",header=FALSE,sep="\t")
mc_CTCF<-mc[rownames(mc)%in%paste0(CTCF_coords$V1,"_",CTCF_coords$V2,"_",CTCF_coords$V3),]
hmc_CTCF<-hmc[rownames(hmc)%in%paste0(CTCF_coords$V1,"_",CTCF_coords$V2,"_",CTCF_coords$V3),]
write.table(mc_CTCF,file="df_5mc_CpG_methyl_CTCF.txt",col.names=TRUE,row.names=TRUE,quote=FALSE,sep="\t")
write.table(hmc_CTCF,file="df_5hmc_CpG_methyl_CTCF.txt",col.names=TRUE,row.names=TRUE,quote=FALSE,sep="\t")

TF_coords<-read.delim("/proj/snic2021-22-81/nobackup/CpG_func_reg_annotated/df_TF_coords.txt",header=FALSE,sep="\t")
mc_TF<-mc[rownames(mc)%in%paste0(TF_coords$V1,"_",TF_coords$V2,"_",TF_coords$V3),]
hmc_TF<-hmc[rownames(hmc)%in%paste0(TF_coords$V1,"_",TF_coords$V2,"_",TF_coords$V3),]
write.table(mc_TF,file="df_5mc_CpG_methyl_TF.txt",col.names=TRUE,row.names=TRUE,quote=FALSE,sep="\t")
write.table(hmc_TF,file="df_5hmc_CpG_methyl_TF.txt",col.names=TRUE,row.names=TRUE,quote=FALSE,sep="\t")

open_chrom_coords<-read.delim("/proj/snic2021-22-81/nobackup/CpG_func_reg_annotated/df_open_chrom_coords.txt",header=FALSE,sep="\t")
mc_open_chrom<-mc[rownames(mc)%in%paste0(open_chrom_coords$V1,"_",open_chrom_coords$V2,"_",open_chrom_coords$V3),]
hmc_open_chrom<-hmc[rownames(hmc)%in%paste0(open_chrom_coords$V1,"_",open_chrom_coords$V2,"_",open_chrom_coords$V3),]
write.table(mc_open_chrom,file="df_5mc_CpG_methyl_open_chrom.txt",col.names=TRUE,row.names=TRUE,quote=FALSE,sep="\t")
write.table(hmc_open_chrom,file="df_5hmc_CpG_methyl_open_chrom.txt",col.names=TRUE,row.names=TRUE,quote=FALSE,sep="\t")



