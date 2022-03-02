# iso_seq_shen

This repository is to cuantify RNA-seq data: </br>
'/home/hawkins/Documents/Cesar/RNA/globus/lordec_reports/rna_seq/rna_seq_fastq' </br>

Data used as reference: </br>
genome1=/home/hawkins/Documents/Cesar/NGSEP/ngsep_tutorial/ZhongmuNo1/ZhongmuNo.1_genome.fasta </br>
transc1=/home/hawkins/Documents/Cesar/NGSEP/ngsep_tutorial/ZhongmuNo1/ZhongmuNo.1.gff3 </br>
Data generated by iso-seq: </br>
genome1='/home/hawkins/Documents/Cesar/NGSEP/ngsep_tutorial/ZhongmuNo1/ZhongmuNo.1_genome.fasta' </br>
transc2='/home/hawkins/Documents/Cesar/RNA/globus/lordec_reports/lordec_trim/bed_Shen/ORF_NMD/blast_corrected_shen.bed' </br>
transc3='/home/hawkins/Documents/Cesar/RNA/globus/lordec_reports/lordec_trim/bed_Shen/ORF_NMD/blast_corrected_shen.gtf' </br>
transc4='/home/hawkins/Documents/Cesar/RNA/globus/lordec_reports/lordec_trim/bed_Shen/ORF_NMD/blast_corrected_shen.2.GFF3' </br>

transc2 has 996554 lines </br>
transc3 has 12970020 lines </br>
transc4 has 9052682 lines  </br>

Run the scripts in the order

1_fastp.sh
This script make a quality control check similar to trimmomatic

2_salmon_quant.sh
This script will aling and quantify your fastq files to obtain a file with tpm in quant.sf

3_mv_salmon.sh
This script moves all quant.sf files in a new directory to import them with tximport R package

tximport.R
