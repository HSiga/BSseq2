setwd("/home/nikolay/WABI/Misc/Zaher/ExprMatrix")
expr<-read.delim("featureCounts_count_OnlyExpressedGenes_TMM_Normalized_coords_and_geneid.txt",header=TRUE,sep="\t")
ensembl_id<-matrix(unlist(strsplit(rownames(expr),"_")),ncol=4,byrow=TRUE)[,1]


library("biomaRt")
ensembl <- useMart("ensembl", dataset="mmusculus_gene_ensembl",host="http://aug2020.archive.ensembl.org")
gene <- getBM(attributes=c('chromosome_name', 'start_position', 'end_position',
                           'ensembl_gene_id','external_gene_name','strand','transcript_start','transcript_end'), 
              filters = 'ensembl_gene_id', values = ensembl_id, mart = ensembl)
head(gene)
write.table(gene,file="ENSEMBL2GENESYMBOL_COORDS_GRCm38.txt",col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")




expr<-read.delim("featureCounts_count_OnlyExpressedGenes_TMM_Normalized_coords_and_geneid.txt",header=TRUE,sep="\t",check.names=FALSE)
head(expr)
my_coords<-matrix(unlist(strsplit(rownames(expr),"_")),ncol=4,byrow=TRUE)
head(my_coords)


ensembl2genesymbol_coords<-read.delim("/home/nikolay/WABI/Misc/Zaher/ExprMatrix/ENSEMBL2GENESYMBOL_COORDS_GRCm38.txt",header=TRUE,sep="\t")
head(ensembl2genesymbol_coords)

intersect_genes<-intersect(my_coords[,1],as.character(ensembl2genesymbol_coords$ensembl_gene_id))
length(intersect_genes)

my_coords<-my_coords[match(intersect_genes,my_coords[,1]),]
dim(my_coords)
head(my_coords)
ensembl2genesymbol_coords<-ensembl2genesymbol_coords[match(intersect_genes,as.character(ensembl2genesymbol_coords$ensembl_gene_id)),]
dim(ensembl2genesymbol_coords)
head(ensembl2genesymbol_coords)

rownames(expr)<-matrix(unlist(strsplit(rownames(expr),"_")),ncol=4,byrow=TRUE)[,1]
head(expr)
dim(expr)
expr<-expr[match(intersect_genes,rownames(expr)),]
dim(expr)
head(expr)
rownames(expr)<-paste0("chr",ensembl2genesymbol_coords$chromosome_name,"_",ensembl2genesymbol_coords$start_position,"_",ensembl2genesymbol_coords$end_position)
head(expr)
gene_length<-as.numeric(ensembl2genesymbol_coords$end_position)-as.numeric(ensembl2genesymbol_coords$start_position)
#hist(gene_length,breaks=100)
min(gene_length)
write.table(expr,file="featureCounts_count_OnlyExpressedGenes_TMM_Normalized_new_coords_GRCm38.txt",col.names = TRUE,row.names = TRUE,quote = FALSE,sep="\t")






pdf("/home/nikolay/WABI/Misc/Zaher/plots/heatmaps_zaher/heatmap_CpG.pdf",paper="a4r",width=297,height=210)
#pdf("/home/nikolay/WABI/Misc/Zaher/plots/heatmaps_zaher/heatmap_dmr.pdf",paper="a4r",width=297,height=210)
#pdf("/home/nikolay/WABI/Misc/Zaher/plots/heatmaps_zaher/heatmap_random1000.pdf", paper="a4r",width=297,height=210)
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
    write.table(aligned_expr,file="/home/nikolay/WABI/Misc/Zaher/plots/heatmaps_zaher/data_points/expr_5mc_intron.txt",col.names=TRUE,row.names=TRUE,quote=FALSE,sep="\t")
    write.table(aligned_methyl,file="/home/nikolay/WABI/Misc/Zaher/plots/heatmaps_zaher/data_points/meth_5mc_intron.txt",col.names=TRUE,row.names=TRUE,quote=FALSE,sep="\t")
    
    coordinates<-as.data.frame(matrix(unlist(strsplit(rownames(aligned_expr),"_")),ncol=6,byrow=TRUE))
    colnames(coordinates)<-c("MethChr","MethStart","MethEnd","GeneChr","GeneStart","GeneEnd")
    head(coordinates)
    #sum(as.numeric(as.character(coordinates$GeneStart))-as.numeric(as.character(coordinates$MethEnd))>0)#promoter
    #sum(as.numeric(as.character(coordinates$MethEnd)) > as.numeric(as.character(coordinates$GeneStart)))#promoter
    sum(as.numeric(as.character(coordinates$GeneStart))-as.numeric(as.character(coordinates$MethStart))>0)#intron
    sum(as.numeric(as.character(coordinates$GeneStart))-as.numeric(as.character(coordinates$MethStart))<0)#intron
    
    gene_length<-as.numeric(matrix(unlist(strsplit(rownames(aligned_expr),"_")),ncol=6,byrow=TRUE)[,6])-
      as.numeric(matrix(unlist(strsplit(rownames(aligned_expr),"_")),ncol=6,byrow=TRUE)[,5])
    min(gene_length)
    max(gene_length)
    #aligned_expr<-aligned_expr[gene_length>50000,]
    #aligned_methyl<-aligned_methyl[gene_length>50000,]
    
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
dev.off()
rho_df<-cbind(rho_df_5mc,rho_df_5hmc)
rho_df
library("pheatmap")
pheatmap(rho_df, display_numbers=TRUE, fontsize=12, main="Spearman correlation rho",cluster_cols=FALSE,cluster_rows=FALSE)
