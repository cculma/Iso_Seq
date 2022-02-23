#!/bin/bash
# commands used in RSEM
source /home/hawkins/miniconda2/bin/activate
conda env list
conda activate alf01

# rsem-prepare-reference rsem1
genome1=/home/hawkins/Documents/Cesar/NGSEP/ngsep_tutorial/ZhongmuNo1/ZhongmuNo.1_genome.fasta
transc1=/home/hawkins/Documents/Cesar/NGSEP/ngsep_tutorial/ZhongmuNo1/ZhongmuNo.1.gff3
rsem-prepare-reference --gff3 $transc1 $genome1 ZhongmuNo.1

mkdir rsem1
mv ZhongmuNo.1.gff3 ZhongmuNo.1.gtf ZhongmuNo.1.grp ZhongmuNo.1.ti ZhongmuNo.1.chrlist ZhongmuNo.1.transcripts.fa ZhongmuNo.1.seq ZhongmuNo.1.idx.fa ZhongmuNo.1.n2g.idx.fa rsem1
cp ZhongmuNo.1_genome.fasta rsem1

# rsem-prepare-reference rsem2
mkdir rsem2
genome1=/home/hawkins/Documents/Cesar/NGSEP/ngsep_tutorial/ZhongmuNo1/ZhongmuNo.1_genome.fasta
transc3=/home/hawkins/Documents/Cesar/RNA/globus/lordec_reports/lordec_trim/bed_Shen/ORF_NMD/blast_corrected_shen.gtf
rsem-prepare-reference --gtf $transc3 $genome1 ZhongmuNo.1
mv ZhongmuNo.1.grp ZhongmuNo.1.ti ZhongmuNo.1.chrlist ZhongmuNo.1.transcripts.fa ZhongmuNo.1.seq ZhongmuNo.1.idx.fa ZhongmuNo.1.n2g.idx.fa rsem2
cp ZhongmuNo.1_genome.fasta rsem2

# STAR-prepare-reference
REF=/home/hawkins/Documents/Cesar/NGSEP/ngsep_tutorial/ZhongmuNo1/ZhongmuNo.1_genome.fasta
DIR=/home/hawkins/Documents/Cesar/NGSEP/ngsep_tutorial/ZhongmuNo1/STAR_index
transc3=/home/hawkins/Documents/Cesar/RNA/globus/lordec_reports/lordec_trim/bed_Shen/ORF_NMD/blast_corrected_shen.gtf
mkdir STAR_index
STAR --runMode genomeGenerate --genomeDir $DIR --genomeFastaFiles $REF --sjdbGTFfile $transc3 --sjdbOverhang 100 

# ReadsAligner using NGSEP
NGSEP=/home/hawkins/Documents/Cesar/NGSEP/ngsep_tutorial/NGSEPcore_4.0.1/NGSEPcore_4.0.1.jar
GENOME=/home/hawkins/Documents/Cesar/NGSEP/ngsep_tutorial/ZhongmuNo1/ZhongmuNo.1_genome.fasta
FM=/home/hawkins/Documents/Cesar/NGSEP/ngsep_tutorial/ZhongmuNo1/ZhongmuNo.1_genome_indexer.fasta
STR=/home/hawkins/Documents/Cesar/NGSEP/ngsep_tutorial/ZhongmuNo1/ZhongmuNo.1_genome_strs.list
java -Xmx280000m -jar ${NGSEP} ReadsAligner -i ${p}.fastq.gz -i2 -o ${p}.bam -r ${GENOME} -d ${FM} -s ${s} -p ILLUMINA -knownSTRs ${STR} -t 40 >& ${p}_ReadsAligner.log