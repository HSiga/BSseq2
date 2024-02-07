# BSseq2
Epigenetics analysis of whole genome Bisulfite (BS) and Oxidative bisulfite (OxBS) sequencing data

This code was developed for the data analysis included in the published paper _¨Epigenetic modulators link mitochondrial redox homeostasis to cardiac function in a sex-dependent manner ¨._

# Brief description of the scripts used for the manuscript:

For data reproducibility, it is recommended to test a down-sampled dataset on a local UNIX-based system (such as MacOS or Linux). For analyzing a large dataset with single-base pair resolution of BS and OxBS sequencing, a high-performance computing (HPC) environment is necessary. Please refer to the specific software/package recommendations for the computational requirements at each step.

The analysis involves the following main steps:

## 1- Pre-process:
The data (BS and oxBS reads) are first cut with Cutadapt (v1.11) and then mapped to the reference genome, here, the mouse genome (GRCm38.p4) using Bismark (v0.19.0). SAMtools (v1.8) can then be used for sorting and indexing. The methylation counts can then be extracted using “bismark_methylation_extractor” tool. 

## 2- Post process:
Post-processing steps involve running different R/bash scripts as follows: 

#### 0_ncpg.R
This code is to calculate the number of CPGs with step size from all CPGs bed files of a reference genome.

#### 1_Pipeline_from_bismarkoutputs.R
The final outputs from previous steps can be used as inputs to this script. From the pre-processing, CytosineReports "...CpG_report.txt.gz", Bismarkcoverage "...bismark.cov.gz" or after being converted to Methylkit format  "..methcounts.gz" can be used as inputs. The script can be run in R (v.3.6) with the main packages. methylKit package (v1.12.0) and MLML2R package (v0.3.3). 

#### 2_5mc-analysis.R
The output from the previous step *cov.gz files are used as input in this step, the scripts can be run in R (v3.4.3) and the main package to use is “bsseq” (v1.14.0 ) 

#### 3_mouse_annotation.bash
Besh script to annotate the DMR files from the previous step.

#### 4_gtf2introns.R:
The script is to extract intron information from mouse GTF-annotation
No main packages are needed

#### 5_PCA_mc_hmc.R:
Building 5mc and 5hmc matrices and making PCA analysis
No main packages are needed

#### 6_IoU.R:
Intersection over union analysis of DMRs for functional elements

#### 7_methyl_per_func_group.R:
Building methylation levels for each group of functional elements

#### 8_select_diff_meth_CpG.R:
Performing differential methylation analysis

#### 9_heatmap.R:
Computing heatmap combining gene expression and DNA methylation information

#### 10_heatmap_pathway.R:
Computing heatmap combining RNAseq pathway analysis with DNA methylation

#### 11_figures.R:
Making figures for the manuscript

