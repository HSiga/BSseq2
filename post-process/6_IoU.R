################################### ANNOTATION ###############################################
cd /proj/snic2021-23-104/bsseq-mlp/202101_analysis/5mc_anaylsis/annotation
sort -k1,1 -k4,4n refGene_mm10.gtf > /proj/snic2021-22-81/nobackup/refGene_mm10.sorted.gtf
bedtools closest -a bsseq_output_dmrs_sorted.bed -b /proj/snic2021-22-81/nobackup/refGene_mm10.sorted.gtf -d > /proj/snic2021-22-81/nobackup/bsseq_output_dmrs_annotated.txt
cd /proj/snic2021-22-81/nobackup
awk '{ if($11 == 0) { print }}' bsseq_output_dmrs_annotated.txt > bsseq_output_dmrs_annotated_d0.txt
awk '{ if($NF == 0) { print }}' mc_hmc_sites_annotated.txt > mc_hmc_sites_annotated_d0.txt



########################## COMPUTING INTERSECTION OVER UNION ##################################
df<-read.delim("5mC_output_dmrs_3WT_vs_3HOM_annotated_d0.txt",header=FALSE,sep="\t")
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



jaccard_5mc<-c(0.0122,0.0757,0.5423,0.0113,0.0061,0.0088,0.0479,0.0255,0.0152,0.0152,0.0109)
jaccard_5mc<-jaccard_5mc/26052
jaccard_5mc
[1] 4.682942e-07 2.905727e-06 2.081606e-05 4.337479e-07 2.341471e-07 3.377860e-07 1.838630e-06 9.788116e-07 5.834485e-07 5.834485e-07 4.183940e-07
 
jaccard_5hmc<-c(0.0604,0.4109,3.0789,0.0505,0.0218,0.0352,0.2386,0.1253,0.0714,0.0729,0.0545)
jaccard_5hmc<-jaccard_5hmc/138369
jaccard_5hmc
[1] 4.365140e-07 2.969596e-06 2.225137e-05 3.649661e-07 1.575497e-07 2.543922e-07 1.724375e-06 9.055497e-07 5.160115e-07 5.268521e-07 3.938744e-07

my_barplot<-rbind(jaccard_5mc,jaccard_5hmc)
barplot(my_barplot,beside=T)


