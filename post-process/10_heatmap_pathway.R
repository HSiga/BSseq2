pathway_zaher<-read.delim("/home/nikolay/WABI/Misc/Zaher/Pathway_IPA_Zaher.txt",header=TRUE,sep="\t")
head(pathway_zaher,2)

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
rho_df<-matrix(NA,nrow=dim(pathway_zaher)[1],ncol=length(func_elements))
for(t in 1:length(func_elements))
{
print(paste0("WORKING WITH ",func_elements[t]," ***************************************************"))

methyl<-read.delim(paste0("/home/nikolay/WABI/Misc/Zaher/",k,"mC_methyl_regions_dmr/df_5",k,"mc_dmr_methyl_",func_elements[t],".txt"),header=TRUE,row.names=1,sep="\t")
methyl<-na.omit(methyl)
colnames(methyl)<-gsub(paste0(k,"mC_"),"",colnames(methyl))
head(methyl)

chr_methyl<-matrix(unlist(strsplit(rownames(methyl),"_")),byrow=TRUE,ncol=3)[,1]
start_methyl<-matrix(unlist(strsplit(rownames(methyl),"_")),byrow=TRUE,ncol=3)[,2]
end_methyl<-matrix(unlist(strsplit(rownames(methyl),"_")),byrow=TRUE,ncol=3)[,3]

expr<-read.delim("/home/nikolay/WABI/Misc/Zaher/ExprMatrix/featureCounts_count_OnlyExpressedGenes_TMM_Normalized_coords_and_genesymbol.txt",
                 header=TRUE,row.names=1,check.names=FALSE,sep="\t")
colnames(expr)<-gsub("MLP-","",colnames(expr))
head(expr)


rho_pathway<-vector()
a<-1:dim(pathway_zaher)[1]; b<-a[seq(0,length(a),10)]
for(j in 1:dim(pathway_zaher)[1])#loop over pathways
{
  pathway_genes<-strsplit(as.character(pathway_zaher$Molecules[j]),split=',', fixed=TRUE)[[1]]
  pathway_genes<-tolower(pathway_genes)
  pathway_genes<-paste(toupper(substr(pathway_genes, 1, 1)), substr(pathway_genes, 2, nchar(pathway_genes)), sep="")
  pathway_genes
  
  expr_pathway<-expr[grepl(paste(pathway_genes,collapse="|"),rownames(expr)),]
  expr_pathway<-expr_pathway[grepl("chrM",rownames(expr_pathway))==FALSE,]
  expr_pathway<-expr_pathway[grepl("chrX",rownames(expr_pathway))==FALSE,]
  head(expr_pathway)
  chr_expr<-matrix(unlist(strsplit(rownames(expr_pathway),"_")),byrow=TRUE,ncol=4)[,2]
  start_expr<-matrix(unlist(strsplit(rownames(expr_pathway),"_")),byrow=TRUE,ncol=4)[,3]
  end_expr<-matrix(unlist(strsplit(rownames(expr_pathway),"_")),byrow=TRUE,ncol=4)[,4]
  
  #Here, for each gene from the pathway we find a matching methylation site
  methyl_pathway<-matrix(NA,ncol=dim(methyl)[2],nrow=dim(expr_pathway)[1])
  methyl_pathway_rownames<-vector()
  for(i in 1:dim(expr_pathway)[1])#loop over genes in each pathway
  {
    methyl_subset<-methyl[chr_methyl==chr_expr[i],]
    start_methyl_subset<-matrix(unlist(strsplit(rownames(methyl_subset),"_")),byrow=TRUE,ncol=3)[,2]
    end_methyl_subset<-matrix(unlist(strsplit(rownames(methyl_subset),"_")),byrow=TRUE,ncol=3)[,3]
    
    methyl_temp<-methyl_subset[which(abs(as.numeric(start_methyl_subset)-as.numeric(start_expr[i]))==
                                       min(abs(as.numeric(start_methyl_subset)-as.numeric(start_expr[i])))),]
    methyl_temp<-methyl_temp[1,]
    methyl_pathway[i,]<-as.numeric(methyl_temp)
    methyl_pathway_rownames<-append(methyl_pathway_rownames,rownames(methyl_temp))
  }
  colnames(methyl_pathway)<-colnames(methyl)
  rownames(methyl_pathway)<-methyl_pathway_rownames
  
  #compute average spearman rho across the 6 samples
  rho<-vector()
  for(s in 1:dim(methyl_pathway)[2])
  {
    if(dim(methyl_pathway)[1]==1)
    {
      rho<-append(rho,NA)
      next
    }
    if(is.na(as.numeric(cor.test(expr_pathway[,s],methyl_pathway[,s],method="spearman")$estimate)))
    {
      rho<-append(rho,NA)
    }
    else
    {
      rho<-append(rho,as.numeric(cor.test(expr_pathway[,s],methyl_pathway[,s],method="spearman")$estimate))
    }
  }
  rho_pathway<-append(rho_pathway,mean(rho,na.rm=TRUE))
  if(j%in%b) print(paste("Finished: ",b[match(j,b)]," pathways",sep=""))
}
rho_pathway[is.na(rho_pathway)]<-0
rho_df[,t]<-rho_pathway
}
colnames(rho_df)<-paste0("5",k,"mc_",func_elements)
rownames(rho_df)<-pathway_zaher$Ingenuity_Canonical_Pathways
if(k=="")
{
  rho_df_5mc<-rho_df
}else
{
  rho_df_5hmc<-rho_df
}
}
rho_df<-cbind(rho_df_5mc,rho_df_5hmc)
rho_df
library("pheatmap")
pheatmap(rho_df, display_numbers=TRUE, fontsize=10, main="Spearman correlation rho",cluster_cols=FALSE,cluster_rows=FALSE)

