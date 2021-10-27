# BSseq2
Analysing whole genome BS and hBS sequencing data

# Brief description of the scripts used for manuscript:
0_ncpg.R: calculate number of CPGs with step size from all CPGs bed file of a genome.

1_Pipeline_from_bismarkoutputs.R: bismark alignment and methylation calling

2_5mc-analysis.R: DMR analysis

3_mouse_annotation.bash: DMR annotation

4_gtf2introns.R: extracting intron infrmation from mouse GTF-annotation

5_PCA_mc_hmc.R: building mc and hmc matrices and making PCA analysis

6_IoU.R: intersection over union analysis of DMRs with respect to functional elements

7_methyl_per_func_group.R: building methylation levels for each group of functional elements

8_select_diff_meth_CpG.R: performing differential methylation analysis

9_heatmap.R: computing heatmap combining gene expression and DNA methylation information

10_heatmap_pathway.R: computing heatmap combining RNAseq pathway analysis with DNA methylation

11_figures.R: making figures for manuscript
