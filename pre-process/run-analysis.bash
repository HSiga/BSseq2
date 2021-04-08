#!/bin/bash


#replace *** with files' initials

for i in ***
do

	bash pipeline_cgex.bash ${i}_R1_merged.fastq.gz ${i}_R2_merged.fastq.gz

done