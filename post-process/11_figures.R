#setwd("/home/nikolay/WABI/Misc/Zaher/")
setwd("/proj/snic2021-22-81/nobackup")

############################## Figure4B, PCA of 5mc and 5hmc ##########################################
library("data.table")
mC_matrix<-fread("mC_matrix.txt",sep="\t")
rownames(mC_matrix)<-mC_matrix$V1
mC_matrix$V1<-NULL
head(mC_matrix)
PC_mC<-prcomp(t(log10(na.omit(mC_matrix) + 1)), center=TRUE, scale=FALSE)
my_color<-ifelse(grepl("HOM",colnames(mC_matrix)),"red","blue")
plot(PC_mC$x[,1:2],main="5mC", xlab="PC1",ylab="PC2",pch=19,col=my_color)
legend("topright",c("HOM","WT"),inset=0.02,fill=c("red","blue"),cex=1.2)
write.table(PC_mC$x,file="Figure4B_5mc.txt",col.names=TRUE,row.names=TRUE,quote=FALSE,sep="\t")

hmC_matrix<-fread("hmC_matrix.txt",sep="\t")
rownames(hmC_matrix)<-hmC_matrix$V1
hmC_matrix$V1<-NULL
head(hmC_matrix)
PC_hmC<-prcomp(t(log10(na.omit(hmC_matrix) + 1)), center=TRUE, scale=FALSE)
my_color<-ifelse(grepl("HOM",colnames(hmC_matrix)),"red","blue")
plot(PC_hmC$x[,1:2],main="5hmC", xlab="PC1",ylab="PC2",pch=19,col=my_color)
legend("bottomright",c("HOM","WT"),inset=0.02,fill=c("red","blue"),cex=1.2)
write.table(PC_hmC$x,file="Figure4B_5hmc.txt",col.names=TRUE,row.names=TRUE,quote=FALSE,sep="\t")


################ Figures 4D&4E the distribution of the 5mC and 5hmC into different regions #############
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
my_barplot<-rbind(mean_HOM_mc_reg,mean_WT_mc_reg)
barplot(my_barplot,beside=T,col=c("blue","orange"),ylab="Methylation Fraction",names=annot,ylim=c(0,1),main="5mC")
legend("topright",c("HOM","WT"),inset=0.02,fill=c("blue","orange"),cex=1.2)
mc_df<-data.frame(mean_HOM_mc_reg=mean_HOM_mc_reg,sd_HOM_mc_reg=sd_HOM_mc_reg,mean_WT_mc_reg=mean_WT_mc_reg,sd_WT_mc_reg=sd_WT_mc_reg)
rownames(mc_df)<-annot
mc_df
write.table(mc_df,file="/home/nikolay/WABI/Misc/Zaher/Manuscript/Figures/Figure4DE/Figure4D.txt",col.names=TRUE,row.names=TRUE,quote=FALSE,sep="\t")

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
my_barplot<-rbind(mean_HOM_hmc_reg,mean_WT_hmc_reg)
barplot(my_barplot,beside=T,col=c("blue","orange"),ylab="Methylation Fraction",names=annot,ylim=c(0,0.1),main="5hmC")
legend("topright",c("HOM","WT"),inset=0.02,fill=c("blue","orange"),cex=1.2)
hmc_df<-data.frame(mean_HOM_hmc_reg=mean_HOM_hmc_reg,sd_HOM_hmc_reg=sd_HOM_hmc_reg,mean_WT_hmc_reg=mean_WT_hmc_reg,sd_WT_hmc_reg=sd_WT_hmc_reg)
rownames(hmc_df)<-annot
hmc_df
write.table(hmc_df,file="/home/nikolay/WABI/Misc/Zaher/Manuscript/Figures/Figure4DE/Figure4E.txt",col.names=TRUE,row.names=TRUE,quote=FALSE,sep="\t")


################ Figures 4F heatmap of correlation between methylation and gene expression #############
methyl_type<-c("","h")
for(k in methyl_type)
{
  if(k=="")
  {
    print("WORKING WITH MC DATA #############################################################################################")
  }
  else
  {
    print("WORKING WITH HMC DATA #############################################################################################")
  }
  func_elements<-c("CDS","exon","intron","promoter","enhancer","CTCF","TF","open_chrom")
  rho_df<-matrix(NA,ncol=3,nrow=length(func_elements))
  for(j in 1:length(func_elements))
  {
    print(paste0("WORKING WITH ",func_elements[j]," ***************************************************"))
    #methyl<-read.delim(paste0("/home/nikolay/WABI/Misc/Zaher/",k,"mC_methyl_regions_diff_meth_CpG/df_5",k,"mc_CpG_methyl_",func_elements[j],".txt"),
    #                   header=TRUE,row.names=1,sep="\t")
    methyl<-read.delim(paste0("/home/nikolay/WABI/Misc/Zaher/",k,"mC_methyl_regions/df_5",k,"mc_methyl_",func_elements[j],".txt"),header=TRUE,row.names=1,sep="\t")
    methyl<-na.omit(methyl)
    colnames(methyl)<-gsub(paste0(k,"mC_"),"",colnames(methyl))
    if(dim(methyl)[1]>5000){methyl<-methyl[sample(1:dim(methyl)[1],5000),]}
    head(methyl)
    
    expr<-read.delim("/home/nikolay/WABI/Misc/Zaher/ExprMatrix/featureCounts_count_OnlyExpressedGenes_TMM_Normalized_new_coords_GRCm38.txt",
                     header=TRUE,row.names=1,check.names=FALSE,sep="\t")
    colnames(expr)<-gsub("MLP-","",colnames(expr))
    head(expr)
    
    
    gene_length<-as.numeric(matrix(unlist(strsplit(rownames(expr),"_")),ncol=3,byrow=TRUE)[,3])-as.numeric(matrix(unlist(strsplit(rownames(expr),"_")),ncol=3,byrow=TRUE)[,2])
    
    #Exclude short genes
    #expr<-expr[gene_length>=10000,]
    
    #Normalization by gene length
    #new_expr<-matrix(NA,ncol=dim(expr)[2],nrow=dim(expr)[1])
    #for(s in 1:length(gene_length))
    #{
    #  new_expr[s,]<-as.numeric(expr[s,])/gene_length[s]
    #}
    #colnames(new_expr)<-colnames(expr)
    #rownames(new_expr)<-rownames(expr)
    #new_expr<-as.data.frame(new_expr)
    #head(new_expr)
    #expr<-new_expr
    
    chr_methyl<-matrix(unlist(strsplit(rownames(methyl),"_")),byrow=TRUE,ncol=3)[,1]
    start_methyl<-matrix(unlist(strsplit(rownames(methyl),"_")),byrow=TRUE,ncol=3)[,2]
    end_methyl<-matrix(unlist(strsplit(rownames(methyl),"_")),byrow=TRUE,ncol=3)[,3]
    
    chr_expr<-matrix(unlist(strsplit(rownames(expr),"_")),byrow=TRUE,ncol=3)[,1]
    start_expr<-matrix(unlist(strsplit(rownames(expr),"_")),byrow=TRUE,ncol=3)[,2]
    end_expr<-matrix(unlist(strsplit(rownames(expr),"_")),byrow=TRUE,ncol=3)[,3]
    
    aligned_expr<-matrix(NA,ncol=dim(expr)[2],nrow=dim(methyl)[1])
    aligned_methyl<-matrix(NA,ncol=dim(methyl)[2],nrow=dim(methyl)[1])
    aligned_rownames<-vector()
    a<-1:dim(methyl)[1]; b<-a[seq(0,length(a),100)]
    for(i in 1:dim(methyl)[1])
    {
      expr_subset<-expr[chr_expr==chr_methyl[i],]
      start_expr_subset<-as.numeric(matrix(unlist(strsplit(rownames(expr_subset),"_")),byrow=TRUE,ncol=3)[,2])
      end_expr_subset<-as.numeric(matrix(unlist(strsplit(rownames(expr_subset),"_")),byrow=TRUE,ncol=3)[,3])
      if(func_elements[j]=="promoter")
      {
        expr_subset<-expr_subset[as.numeric(start_expr_subset)-as.numeric(end_methyl[i])>0,]
      }
      if(func_elements[j]=="intron")
      {
        expr_subset<-expr_subset[as.numeric(start_expr_subset)-as.numeric(start_methyl[i])<0,]
      }
      if(dim(expr_subset)[1]==0)
      {
        aligned_rownames<-append(aligned_rownames,NA)
        next
      }
      start_expr_subset<-as.numeric(matrix(unlist(strsplit(rownames(expr_subset),"_")),byrow=TRUE,ncol=3)[,2])
      end_expr_subset<-as.numeric(matrix(unlist(strsplit(rownames(expr_subset),"_")),byrow=TRUE,ncol=3)[,3])
      expr_temp<-expr_subset[which(abs(as.numeric(start_expr_subset)-as.numeric(end_methyl[i]))==min(abs(as.numeric(start_expr_subset)-as.numeric(end_methyl[i])))),]
      
      expr_temp<-expr_temp[1,]
      expr_reg<-as.numeric(expr_temp)
      expr_reg_HOM<-as.numeric(subset(expr_temp,select=colnames(expr_temp)[grepl("HOM",colnames(expr_temp))]))
      expr_reg_WT<-as.numeric(subset(expr_temp,select=colnames(expr_temp)[grepl("WT",colnames(expr_temp))]))
      head(expr_temp)
      
      methyl_temp<-methyl[i,]
      methyl_reg<-as.numeric(methyl_temp)
      methyl_reg_HOM<-as.numeric(subset(methyl_temp,select=colnames(methyl_temp)[grepl("HOM",colnames(methyl_temp))]))
      methyl_reg_WT<-as.numeric(subset(methyl_temp,select=colnames(methyl_temp)[grepl("WT",colnames(methyl_temp))]))
      head(methyl_temp)
      
      aligned_rownames<-append(aligned_rownames,paste0(rownames(methyl_temp),"_",rownames(expr_temp)))
      aligned_expr[i,]<-as.numeric(expr_temp)
      aligned_methyl[i,]<-as.numeric(methyl_temp)
      
      if(i%in%b) print(paste("Finished: ",b[match(i,b)]," methylation sites",sep=""))
    }
    colnames(aligned_methyl)<-colnames(methyl)
    colnames(aligned_expr)<-colnames(expr)
    rownames(aligned_methyl)<-aligned_rownames
    rownames(aligned_expr)<-aligned_rownames
    
    head(aligned_expr)
    head(aligned_methyl)
    
    coordinates<-as.data.frame(matrix(unlist(strsplit(rownames(aligned_expr),"_")),ncol=6,byrow=TRUE))
    colnames(coordinates)<-c("MethChr","MethStart","MethEnd","GeneChr","GeneStart","GeneEnd")
    head(coordinates)
    sum(as.numeric(as.character(coordinates$GeneStart))-as.numeric(as.character(coordinates$MethStart))>0)#intron
    sum(as.numeric(as.character(coordinates$GeneStart))-as.numeric(as.character(coordinates$MethStart))<0)#intron
    
    gene_length<-as.numeric(matrix(unlist(strsplit(rownames(aligned_expr),"_")),ncol=6,byrow=TRUE)[,6])-
      as.numeric(matrix(unlist(strsplit(rownames(aligned_expr),"_")),ncol=6,byrow=TRUE)[,5])
    min(gene_length)
    max(gene_length)
    
    par(mfrow=c(2,3),oma = c(0, 0, 2, 0))
    plot(log10(aligned_expr[,1]+1)~log10(aligned_methyl[,1]+1),xlab="log10( Methylation )",ylab="log10( Gene Expression )",main="HOM2355")
    abline(lm(log10(aligned_expr[,1]+1)~log10(aligned_methyl[,1]+1)))
    mtext(paste0("Spearman rho = ",as.numeric(cor.test(aligned_expr[,1],aligned_methyl[,1],method="spearman")$estimate)))
    plot(log10(aligned_expr[,2]+1)~log10(aligned_methyl[,2]+1),xlab="log10( Methylation )",ylab="log10( Gene Expression )",main="HOM2356")
    abline(lm(log10(aligned_expr[,2]+1)~log10(aligned_methyl[,2]+1)))
    mtext(paste0("Spearman rho = ",as.numeric(cor.test(aligned_expr[,2],aligned_methyl[,2],method="spearman")$estimate)))
    plot(log10(aligned_expr[,3]+1)~log10(aligned_methyl[,3]+1),xlab="log10( Methylation )",ylab="log10( Gene Expression )",main="HOM2357")
    abline(lm(log10(aligned_expr[,3]+1)~log10(aligned_methyl[,3]+1)))
    mtext(paste0("Spearman rho = ",as.numeric(cor.test(aligned_expr[,3],aligned_methyl[,3],method="spearman")$estimate)))
    plot(log10(aligned_expr[,4]+1)~log10(aligned_methyl[,4]+1),xlab="log10( Methylation )",ylab="log10( Gene Expression )",main="WT2")
    abline(lm(log10(aligned_expr[,4]+1)~log10(aligned_methyl[,4]+1)))
    mtext(paste0("Spearman rho = ",as.numeric(cor.test(aligned_expr[,4],aligned_methyl[,4],method="spearman")$estimate)))
    plot(log10(aligned_expr[,5]+1)~log10(aligned_methyl[,5]+1),xlab="log10( Methylation )",ylab="log10( Gene Expression )",main="WT3")
    abline(lm(log10(aligned_expr[,5]+1)~log10(aligned_methyl[,5]+1)))
    mtext(paste0("Spearman rho = ",as.numeric(cor.test(aligned_expr[,5],aligned_methyl[,5],method="spearman")$estimate)))
    plot(log10(aligned_expr[,6]+1)~log10(aligned_methyl[,6]+1),xlab="log10( Methylation )",ylab="log10( Gene Expression )",main="WT3B")
    abline(lm(log10(aligned_expr[,6]+1)~log10(aligned_methyl[,6]+1)))
    mtext(paste0("Spearman rho = ",as.numeric(cor.test(aligned_expr[,6],aligned_methyl[,6],method="spearman")$estimate)))
    mtext(paste0("5",k,"mc:"," ",func_elements[j]), outer = TRUE, cex = 1.5)
    
    spear<-function(x){return(as.numeric(cor.test(aligned_expr[,x],aligned_methyl[,x],method="spearman")$estimate))}
    
    rho<-c(spear(1),spear(2),spear(3),spear(4),spear(5),spear(6))
    rho_HOM<-c(spear(1),spear(2),spear(3))
    rho_WT<-c(spear(4),spear(5),spear(6))
    rho_df[j,]<-c(mean(rho_WT,na.rm=TRUE), mean(rho_HOM,na.rm=TRUE), mean(rho,na.rm=TRUE))
  }
  colnames(rho_df)<-c(paste0("WT_5",k,"mc"),paste0("HOM_5",k,"mc"),paste0("Combined_5",k,"mc"))
  rownames(rho_df)<-func_elements
  if(k=="")
  {
    rho_df_5mc<-rho_df
  }
  else
  {
    rho_df_5hmc<-rho_df
  }
}
rho_df<-cbind(rho_df_5mc,rho_df_5hmc)
rho_df
write.table(rho_df,file="/home/nikolay/WABI/Misc/Zaher/Manuscript/Figures/Figure4F.txt",col.names=TRUE,row.names=TRUE,quote=FALSE,sep="\t")
library("pheatmap")
pheatmap(rho_df, display_numbers=TRUE, fontsize=12, main="Spearman correlation rho",cluster_cols=FALSE,cluster_rows=FALSE)


################ Figure 4G DMR intersection with functional regions #############
df<-read.delim("/home/nikolay/WABI/Misc/Zaher/DMR_5mc_5hmc/5hmC_output_dmrs_3WT_vs_3HOM.1perc_diff_annotated_d0.txt",header=FALSE,sep="\t")
df_total<-read.delim("/home/nikolay/WABI/Misc/Zaher/reference_data/mm10_annot_nodups.sorted.gtf",header=FALSE,sep="\t")

intersect_exon<-dim(df[as.character(df$V6)=="exon",])[1]
print(intersect_exon)
union_exon<-dim(df_total[as.character(df_total$V3)=="exon",])[1]
print(union_exon)
intersect_exon/union_exon


intersect_intron<-dim(df[as.character(df$V6)=="intron",])[1]
print(intersect_intron)
union_intron<-dim(df_total[as.character(df_total$V3)=="intron",])[1]
print(union_intron)
intersect_intron/union_intron


intersect_gene_body<-dim(df[as.character(df$V6)=="transcript",])[1]
print(intersect_gene_body)
union_gene_body<-dim(df_total[as.character(df_total$V3)=="transcript",])[1]
print(union_gene_body)
intersect_gene_body/union_gene_body


intersect_CDS<-dim(df[as.character(df$V6)=="CDS",])[1]
print(intersect_CDS)
union_CDS<-dim(df_total[as.character(df_total$V3)=="CDS",])[1]
print(union_CDS)
intersect_CDS/union_CDS


intersect_start_codon<-dim(df[as.character(df$V6)=="start_codon",])[1]
print(intersect_start_codon)
union_start_codon<-dim(df_total[as.character(df_total$V3)=="start_codon",])[1]
print(union_start_codon)
intersect_start_codon/union_start_codon


intersect_stop_codon<-dim(df[as.character(df$V6)=="stop_codon",])[1]
print(intersect_stop_codon)
union_stop_codon<-dim(df_total[as.character(df_total$V3)=="stop_codon",])[1]
print(union_stop_codon)
intersect_stop_codon/union_stop_codon


intersect_promoter<-dim(df[grepl("Promoter",as.character(df$V12)),])[1]
print(intersect_promoter)
union_promoter<-dim(df_total[grepl("Promoter",as.character(df_total$V9)),])[1]
print(union_promoter)
intersect_promoter/union_promoter



intersect_enhancer<-dim(df[grepl("Enhancer",as.character(df$V12)),])[1]
print(intersect_enhancer)
union_enhancer<-dim(df_total[grepl("Enhancer",as.character(df_total$V9)),])[1]
print(union_enhancer)
intersect_enhancer/union_enhancer


intersect_CTCF<-dim(df[grepl("CTCF",as.character(df$V12)),])[1]
print(intersect_CTCF)
union_CTCF<-dim(df_total[grepl("CTCF",as.character(df_total$V9)),])[1]
print(union_CTCF)
intersect_CTCF/union_CTCF



intersect_TF<-dim(df[grepl("TF",as.character(df$V12)),])[1]
print(intersect_TF)
union_TF<-dim(df_total[grepl("TF",as.character(df_total$V9)),])[1]
print(union_TF)
intersect_TF/union_TF


intersect_open_chromatin<-dim(df[grepl("Open",as.character(df$V12)),])[1]
print(intersect_open_chromatin)
union_open_chromatin<-dim(df_total[grepl("Open",as.character(df_total$V9)),])[1]
print(union_open_chromatin)
intersect_open_chromatin/union_open_chromatin

#jaccard_5mc<-c(0.0122,0.0757,0.5423,0.0113,0.0061,0.0088,0.0479,0.0255,0.0152,0.0152,0.0109)
jaccard_5mc<-c(intersect_exon/union_exon,
               intersect_intron/union_intron,
               intersect_CDS/union_CDS,
               intersect_start_codon/union_start_codon,
               intersect_stop_codon/union_stop_codon,
               intersect_promoter/union_promoter,
               intersect_enhancer/union_enhancer,
               intersect_CTCF/union_CTCF,
               intersect_TF/union_TF,
               intersect_open_chromatin/union_open_chromatin)
jaccard_5mc<-jaccard_5mc/26052 #total number of 5mc DMRs (10% difference between WT and HOM)
jaccard_5mc

#jaccard_5hmc<-c(0.0604,0.4109,3.0789,0.0505,0.0218,0.0352,0.2386,0.1253,0.0714,0.0729,0.0545)
jaccard_5hmc<-c(intersect_exon/union_exon,
               intersect_intron/union_intron,
               intersect_CDS/union_CDS,
               intersect_start_codon/union_start_codon,
               intersect_stop_codon/union_stop_codon,
               intersect_promoter/union_promoter,
               intersect_enhancer/union_enhancer,
               intersect_CTCF/union_CTCF,
               intersect_TF/union_TF,
               intersect_open_chromatin/union_open_chromatin)
jaccard_5hmc<-jaccard_5hmc/138369 #total number of 5hmc DMRs (1% difference between WT and HOM)
jaccard_5hmc

my_barplot<-rbind(jaccard_5mc,jaccard_5hmc)
colnames(my_barplot)<-c("exon","intron","CDS","start","stop","promoter","enhancer","CTCF","TF","OpenChrom")
my_barplot
par(mfrow=c(1,1))
barplot(my_barplot,beside=T,ylab="Intersection Over Union (IoU)",
        names=c("exon","intron","CDS","start","stop","promoter","enhancer","CTCF","TF","OpenChrom"),col=c("darkgreen","red"),cex.axis=1,
        cex.names=1.5)
legend("topright",c("5MC","5HMC"),inset=0.02,fill=c("darkgreen","red"),cex=2)
write.table(my_barplot,file="/home/nikolay/WABI/Misc/Zaher/Manuscript/Figures/Figure4G.txt",col.names=TRUE,row.names=TRUE,quote=FALSE,sep="\t")

######################## Figure 4H list of genes with DMRs in their intronic regions for 5hmc ###############################
setwd("/home/nikolay/WABI/Misc/Zaher/DMR_5mc_5hmc")
dmr<-read.delim("5hmC_output_dmrs_3WT_vs_3HOM.1perc_diff.txt",header=TRUE,sep="\t")
dmr<-dmr[order(-abs(dmr$meanDiff)),]
head(dmr)

dmr_annot<-read.delim("5hmC_output_dmrs_3WT_vs_3HOM.1perc_diff_annotated_d0.txt",header=FALSE,sep="\t")
dmr_annot<-dmr_annot[as.character(dmr_annot$V6)=="intron",]
head(dmr_annot)

dmr<-dmr[match(paste0(dmr_annot$V1,"_",dmr_annot$V2,"_",dmr_annot$V3),paste0(dmr$chr,"_",dmr$start,"_",dmr$end)),]
dmr<-dmr[order(-abs(dmr$meanDiff)),]
head(dmr)
my_dmr<-dmr[1:5000,]
head(my_dmr)
tail(my_dmr)

my_dmr_annot<-dmr_annot[match(paste0(my_dmr$chr,"_",my_dmr$start,"_",my_dmr$end),paste0(dmr_annot$V1,"_",dmr_annot$V2,"_",dmr_annot$V3)),]
my_dmr_annot$meanDiff<-my_dmr$meanDiff
head(my_dmr_annot)
head(gsub(";","",matrix(unlist(strsplit(as.character(my_dmr_annot$V12)," gene_name ")),ncol=2,byrow=TRUE)[,2]))
my_output<-data.frame(genes=gsub(";","",matrix(unlist(strsplit(as.character(my_dmr_annot$V12)," gene_name ")),ncol=2,byrow=TRUE)[,2]),meanDiff=my_dmr_annot$meanDiff)
head(my_output)
my_output<-my_output[!duplicated(as.character(my_output$genes)),]
head(my_output)
dim(my_output)
write.table(my_output,file="/home/nikolay/WABI/Misc/Zaher/Manuscript/Figures/Figure4H.txt",col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")

