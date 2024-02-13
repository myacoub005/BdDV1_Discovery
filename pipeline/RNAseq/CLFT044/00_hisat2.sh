#!/bin/bash
#SBATCH -p intel --time 12-0:00:00 --ntasks 8 --nodes 1 --mem 24G --out logs/hisat.log

module load hisat2
module load samtools

#bowtie2-build Vir_coverage/Batrachochytrium_dendrobatidis_CLFT044.scaffolds.fa CLFT044

#bowtie2 -x CLFT044 -1 fastq/C44_A_1.fq.gz -2 fastq/C44_A_2.fq.gz | samtools view -bS - > CLFT044.bam

#bowtie2 --bam CLFT044 -1 reads_1.fastq -2 reads_2.fastq output.bam

$OUTDIR=hisat

#hisat2-build Vir_coverage/Batrachochytrium_dendrobatidis_CLFT044.scaffolds.fa CLFT044

#hisat2 -x CLFT044 -1 fastq/C44_A_1.fq.gz -2 fastq/C44_A_2.fq.gz -S CLFT044.bam

samtools sort CLFT044.bam > CLFT044.sorted.bam

samtools mpileup CLFT044.sorted.bam | awk '{print $1"\t"$2"\t"$4}' > $OUTDIR/coverage.res.tsv
