#!/bin/bash
for i in *fastq.gz

do
	p=$(echo $i| cut -d'_' -f 1,2);
	echo ${p};
#	fastp -i ${i} -o out_${i};
	fastp -i ${p}_R1_001.fastq.gz -I ${p}_R2_001.fastq.gz -o out_${p}_R1_001.fastq.gz -O out_${p}_R2_001.fastq.gz;
done
