#!/bin/bash
salmon=/home/hawkins/Documents/Cesar/RNA/Iso_assay/salmon-1.7.0_linux_x86_64/bin/salmon;
salmon_index=/home/hawkins/Documents/Cesar/NGSEP/ngsep_tutorial/ZhongmuNo1/salmon_index/salmon_index/;

for i in *R1_001.fastq.gz;
do
  	p=$(echo $i| cut -d'_' -f 1,2);
	echo "Processing sample ${p}";
	$salmon quant -i $salmon_index -l A -1 ${p}_R1_001.fastq.gz -2 ${p}_R2_001.fastq.gz --validateMappings -o salmon_shen_${p};
done
