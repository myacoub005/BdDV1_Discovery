#!/usr/bin/bash -l
#SBATCH -p short -N 1 -n 32 --mem 32gb --out logs/mapVirus_nanopore.log
CPU=1
if [ $SLURM_CPUS_ON_NODE ]; then
  CPU=$SLURM_CPUS_ON_NODE
fi

module load minimap2
module load samtools
OUT=analysis/map_virus
READS=CLFT067.VirusReads.fq
VIRGENOME=genome/Batrachochytrium_dendrobatidis_CLFT044.Polished.scaffolds.fa
mkdir -p $OUT
minimap2 -t $CPU -x map-ont -a $VIRGENOME $READS > $OUT/CLFT67_to_Virus.sam
minimap2 -t $CPU -x map-ont --cs=long $VIRGENOME $READS > $OUT/CLFT67_to_Virus.paf

samtools sort -O BAM -o $OUT/CLFT67_to_Virus.bam --threads $CPU $OUT/CLFT67_to_Virus.sam
samtools index  $OUT/CLFT67_to_Virus.bam
samtools view -O BAM -o $OUT/CLFT67_to_Virus.maponly.bam -F 4 --threads $CPU $OUT/CLFT67_to_Virus.bam
samtools index $OUT/CLFT67_to_Virus.maponly.bam

# maps virus reads back to reference genome
