#!/bin/bash
#SBATCH -p short -N 1 -n 16 --mem 16gb --out spades.%a.log


module load spades

CPU=$SLURM_CPUS_ON_NODE
N=${SLURM_ARRAY_TASK_ID}

if [ ! $N ]; then
    N=$1
    if [ ! $N ]; then
        echo "Need an array id or cmdline val for the job"
        exit
    fi
fi

INDIR=Virus_fq
OUT=BdDV1_Illumina_assemblies
SAMPLEFILE=samples_PE.dat

if [ ! -f $REF.bwt ]; then
    bwa index $REF
fi

mkdir -p $OUT

sed -n ${N}p $SAMPLEFILE | while read LIBRARY STRAIN
do

spades.py --only-assembler -s $INDIR/${STRAIN}.fq -o $OUT/${STRAIN}

done
