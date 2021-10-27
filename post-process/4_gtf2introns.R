setwd("/home/nikolay/WABI/Misc/Zaher")

df<-read.delim("refGene_mm10.gtf",header=FALSE,sep="\t")
head(df,8)
levels(as.factor(df$V3))

df_CDS<-df[as.character(df$V3)=="CDS",]
df_start_codon<-df[as.character(df$V3)=="start_codon",]
df_stop_codon<-df[as.character(df$V3)=="stop_codon",]
df_exon<-df[as.character(df$V3)=="exon",]
df_gene<-df[as.character(df$V3)=="transcript",]
df_gene_exon<-df[as.character(df$V3)=="transcript" | as.character(df$V3)=="exon",]

n_genes<-sum(as.character(df_gene_exon$V3)=="transcript")
gene_index<-as.numeric(rownames(df_gene_exon)[as.character(df_gene_exon$V3)=="transcript"])
df_intron<-data.frame(matrix(ncol=dim(df_gene_exon)[2]))
colnames(df_intron)<-colnames(df_gene_exon)
df_intron<-df_intron[-1,]
for(i in 1:n_genes)
{
  print(paste0("WORKING WITH GENE ",i))
  gene_start<-df_gene_exon$V4[i]
  gene_end<-df_gene_exon$V5[i]
  df_gene_temp<-df_gene_exon[as.numeric(rownames(df_gene_exon))>gene_index[i] & as.numeric(rownames(df_gene_exon))<gene_index[i+1],]
  
  df_intron_temp<-data.frame(matrix(ncol=dim(df_gene_temp)[2]))
  colnames(df_intron_temp)<-colnames(df_gene_temp)
  df_intron_temp<-df_intron_temp[-1,]
  for(j in 1:dim(df_gene_temp)[1])
  {
    df_intron_temp<-rbind(df_intron_temp,df_gene_temp[j,])
    df_intron_temp[j,]$V4<-df_intron_temp[j,]$V5
    df_intron_temp[j,]$V5<-df_gene_temp[j+1,]$V4
  }
  df_intron_temp<-df_intron_temp[-nrow(df_intron_temp),]
  df_intron_temp<-data.frame(lapply(df_intron_temp,function(x){gsub("exon","intron",x)}))
  df_intron<-rbind(df_intron,df_intron_temp)
}
