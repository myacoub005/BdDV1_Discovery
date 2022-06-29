#!/bin/bash
#SBATCH -p short -N 1 -n 16 --mem 16gb --out logs/mpileup.log

module load samtools

samtools mpileup analysis/map_virus/CLFT67_to_Virus.bam --output bedfiles/CLFT067_to_Virus.bam

awk '{print $1,$2,$4}' bedfiles/CLFT067_to_Virus.bam > CLFT067_to_Virus.bed

# gives locations to where the virus reads overlaps with reference genome
