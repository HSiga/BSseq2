#PREPARE BED-FILE WITH DMR REGIONS
cut -f1-3 ../5mC_output_dmrs_3WT_vs_3HOM.txt > bsseq_output_dmrs.bed
sed 1d bsseq_output_dmrs.bed > bsseq_output_dmrs_temp.bed
rm bsseq_output_dmrs.bed
mv bsseq_output_dmrs_temp.bed bsseq_output_dmrs.bed
 
#DOWNLOAD mm10 GTF ANNOTATION AND PREPROCESS IT (REMOVE "CHR", TRANSCRITS, CDS AND DUPLICATES)
#Download genePredToGTF from here http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/
module load ucsc-utilities

mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -N -e "select * from refGene;" mm10 | cut -f2- | genePredToGtf file stdin refGene.gtf

mv refGene.gtf refGene_mm10.gtf
#Skip removing chr
#sed 's/chr//g' refGene_mm10.gtf > refGene_mm10.gtf

grep -wv transcript refGene_mm10.gtf > refGene_mm10_notranscript.gtf
grep -wv CDS refGene_mm10_notranscript.gtf > refGene_mm10_notranscript_noCDS.gtf
 
#NOW REMOVE DUPLICATE ANNOTATIONS USING R
df<-read.delim("refGene_mm10_notranscript_noCDS.gtf",header=FALSE)
df$coord<-paste0(df$V1,"_",df$V4,"_",df$V5)
df_nodups<-df[!duplicated(df$coord),]
df_nodups$coord<-NULL
write.table(df_nodups,file="refGene_mm10_notranscript_noCDS_nodups.gtf",col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")

#Sort both
bedtools sort -i bsseq_output_dmrs.bed > bsseq_output_dmrs_sorted.bed
bedtools sort -i refGene_mm10_notranscript_noCDS_nodups.gtf > refGene_mm10_notranscript_noCDS_nodups_sorted.gtf

#RUN ANNOTATION WITH BEDTOOLS
bedtools closest -a bsseq_output_dmrs_sorted.bed -b refGene_mm10_notranscript_noCDS_nodups_sorted.gtf -d > bsseq_output_dmrs_annotated.txt
 
#MAKE OUTPUT PRETTY IN R
df<-read.delim("../5mC_output_dmrs_3WT_vs_3HOM.txt",header=TRUE)
annot<-read.delim("bsseq_output_dmrs_annotated.txt",header=FALSE)
annot<-annot[match(paste0(df$chr,"_",df$start,"_",df$end),paste0(annot$V1,"_",annot$V2,"_",annot$V3)),]
df$gene<-annot$V12[match(paste0(df$chr,"_",df$start,"_",df$end),paste0(annot$V1,"_",annot$V2,"_",annot$V3))]
df$dist2gene<-annot$V13[match(paste0(df$chr,"_",df$start,"_",df$end),paste0(annot$V1,"_",annot$V2,"_",annot$V3))]
write.table(df,file="bsseq_output_dmrs_annotated.txt",col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")