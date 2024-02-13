#!/usr/bin/bash
#SBATCH --nodes 1 --ntasks 24 --mem 24G -p short -J bwacov --out logs/bwa.log --time 2:00:00

module load samtools
module load bwa

CPU=$SLURM_CPUS_ON_NODE

INDIR=fastq
BASE=C44_A
REF=Vir_coverage/Batrachochytrium_dendrobatidis_CLFT044.scaffolds.fa
OUTDIR=bwa

bwa index $REF

bwa mem -B 40 -O 60 -E 10 -L 100 -t $CPU $REF $INDIR/${BASE}_1.fq.gz $INDIR/${BASE}_2.fq.gz > $OUTDIR/${BASE}.sam

samtools fixmate -O bam $OUTDIR/${BASE}.sam $OUTDIR/${BASE}.bam

samtools sort $OUTDIR/${BASE}.bam > $OUTDIR/${BASE}.sorted.bam

samtools mpileup $OUTDIR/${BASE}.sorted.bam | awk '{print $1"\t"$2"\t"$4}' > $OUTDIR/coverage.res.tsv

#samtools sort CLFT044.bam > CLFT044.sorted.bam

#samtools mpileup CLFT044.sorted.bam | awk '{print $1"\t"$2"\t"$4}' > coverage.res.tsv
