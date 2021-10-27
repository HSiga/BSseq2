library("data.table")
setwd("/crex/proj/snic2021-23-104/bsseq-mlp/202101_analysis")
samples<-list.files(pattern="cov.txt.gz")
methyl<-list()
for(i in 1:length(samples))
{
  print(paste0("READING SAMPLE ",samples[i]))
  methyl[[i]]<-as.data.frame(fread(paste0("zcat ",samples[i]),header=TRUE))
  rownames(methyl[[i]])<-paste0(methyl[[i]]$chr,"_",methyl[[i]]$start,"_",methyl[[i]]$end)
  methyl[[i]]$chr<-NULL
  methyl[[i]]$start<-NULL
  methyl[[i]]$end<-NULL
  methyl[[i]]$strand<-NULL
}
coord_intersect<-Reduce(intersect,list(rownames(methyl[[1]]),rownames(methyl[[2]]),rownames(methyl[[3]]),rownames(methyl[[4]]),rownames(methyl[[5]]),rownames(methyl[[6]])))

for(i in 1:length(samples))
{
  print(paste0("PRUNING SAMPLE ",samples[i]))
  methyl[[i]]<-methyl[[i]][match(coord_intersect,rownames(methyl[[i]])),]
}
methyl_matrix<-Reduce(cbind,list(methyl[[1]],methyl[[2]],methyl[[3]],methyl[[4]],methyl[[5]],methyl[[6]]))
cov_matrix<-subset(methyl_matrix,select=colnames(methyl_matrix)[grepl("COVbs_oxbs_sample",colnames(methyl_matrix))])
mC_matrix<-subset(methyl_matrix,select=colnames(methyl_matrix)[grepl("^mC",colnames(methyl_matrix))])
hmC_matrix<-subset(methyl_matrix,select=colnames(methyl_matrix)[grepl("hmC",colnames(methyl_matrix))])

write.table(methyl_matrix,file="/proj/snic2021-22-81/nobackup/methyl_matrix.txt",col.names=TRUE,row.names=TRUE,quote=FALSE,sep="\t")
write.table(cov_matrix,file="/proj/snic2021-22-81/nobackup/cov_matrix.txt",col.names=TRUE,row.names=TRUE,quote=FALSE,sep="\t")
write.table(mC_matrix,file="/proj/snic2021-22-81/nobackup/mC_matrix.txt",col.names=TRUE,row.names=TRUE,quote=FALSE,sep="\t")
write.table(hmC_matrix,file="/proj/snic2021-22-81/nobackup/hmC_matrix.txt",col.names=TRUE,row.names=TRUE,quote=FALSE,sep="\t")

pdf("/proj/snic2021-22-81/nobackup/hist_cov.pdf")
hist(log10(as.numeric(rowMeans(cov_matrix))+1),breaks=100,xlab="log10(Mean Coverage)")
dev.off()


PC_mC<-prcomp(t(log10(na.omit(mC_matrix) + 1)), center=TRUE, scale=FALSE)
my_color<-ifelse(grepl("HOM",colnames(mC_matrix)),"red","blue")
pdf("/proj/snic2021-22-81/nobackup/PCA_5mc.pdf")
plot(PC_mC$x[,1:2],main="5mC", xlab="PC1",ylab="PC2",pch=19,col=my_color)
legend("topright",c("HOM","WT"),inset=0.02,fill=c("red","blue"),cex=1.2)
dev.off()


PC_hmC<-prcomp(t(log10(na.omit(hmC_matrix) + 1)), center=TRUE, scale=FALSE)
my_color<-ifelse(grepl("HOM",colnames(hmC_matrix)),"red","blue")
pdf("/proj/snic2021-22-81/nobackup/PCA_5hmc.pdf")
plot(PC_hmC$x[,1:2],main="5hmC", xlab="PC1",ylab="PC2",pch=19,col=my_color)
legend("bottomright",c("HOM","WT"),inset=0.02,fill=c("red","blue"),cex=1.2)
dev.off()

