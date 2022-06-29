#!/bin/bash
#SBATCH -p short -N 1 -n 16 --mem 16gb --out logs/samtools.fq.log

samtools view -h -b analysis/map_longread/CLFT67_to_VirCtg.bam "contig_138:3537-7941" > CLFT067.Virus.bam
samtools fastq CLFT067.Virus.bam > CLFT067.VirusReads.fq

##writes virus reads out as new fq data
