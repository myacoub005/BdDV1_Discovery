#!/bin/bash
#SBATCH -p short -N 1 -n 16 --mem 16gb --out prodigal.%a.log


module load prodigal

CPU=$SLURM_CPUS_ON_NODE
N=${SLURM_ARRAY_TASK_ID}

if [ ! $N ]; then
    N=$1
    if [ ! $N ]; then
        echo "Need an array id or cmdline val for the job"
        exit
    fi
fi

INDIR=asm_unassembled
SAMPLEFILE=samples_PE.dat

if [ ! -f $REF.bwt ]; then
    bwa index $REF
fi


sed -n ${N}p $SAMPLEFILE | while read LIBRARY STRAIN
do

prodigal -i $INDIR/${STRAIN}/scaffolds.fasta -o $INDIR/${STRAIN}/${STRAIN} -d $INDIR/${STRAIN}/${STRAIN}.orfs.fasta -a ${INDIR}/${STRAIN}/${STRAIN}.proteins.fasta

done
