#!/usr/bin/bash -l
#SBATCH -p short -N 1 -n 32 --mem 32gb --out logs/maplongread_nanopore.log
CPU=1
if [ $SLURM_CPUS_ON_NODE ]; then
  CPU=$SLURM_CPUS_ON_NODE
fi

module load minimap2
module load samtools
OUT=analysis/map_longread
READS=/rhome/myaco005/shared/projects/Chytrid/BdBRAZIL_Nanopore/nanopore_fq/67_nanopore.fastq.gz
VIRGENOME=genome/Batrachochytrium_dendrobatidis_CLFT044.Polished.scaffolds.fa
mkdir -p $OUT
minimap2 -t $CPU -x map-ont -a $VIRGENOME $READS > $OUT/CLFT67_to_VirCtg.sam
minimap2 -t $CPU -x map-ont --cs=long $VIRGENOME $READS > $OUT/CLFT67_to_VirCtg.paf

samtools sort -O BAM -o $OUT/CLFT67_to_VirCtg.bam --threads $CPU $OUT/CLFT67_to_VirCtg.sam
samtools index  $OUT/CLFT67_to_VirCtg.bam
samtools view -O BAM -o $OUT/CLFT67_to_VirCtg.maponly.bam -F 4 --threads $CPU $OUT/CLFT67_to_VirCtg.bam
samtools index $OUT/CLFT67_to_VirCtg.maponly.bam

#### Gives you reads mapped only to viral genome from nanopore ONT data ####
