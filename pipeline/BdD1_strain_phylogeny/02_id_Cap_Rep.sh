#!/bin/bash
#SBATCH -p short -N 1 -n 16 --mem 16gb --out logs/tblast.%a.log


module load ncbi-blast

CPU=$SLURM_CPUS_ON_NODE
N=${SLURM_ARRAY_TASK_ID}

if [ ! $N ]; then
    N=$1
    if [ ! $N ]; then
        echo "Need an array id or cmdline val for the job"
        exit
    fi
fi

INDIR=BdDV1_Illumina_assemblies
SAMPLEFILE=samples_PE.dat
OUT=TBLAST_TEST

if [ ! -f $REF.bwt ]; then
    bwa index $REF
fi


sed -n ${N}p $SAMPLEFILE | while read LIBRARY STRAIN
do
mkdir $OUT

tblastn -query BdDV1.REP.fasta -subject $INDIR/${STRAIN}/${STRAIN}.orfs.fasta -evalue 0.00001 -outfmt 6 -max_hsps 1 -max_target_seqs 1 -out $OUT/${STRAIN}.REP.out
tblastn -query BdDV1.CAP.fasta -subject $INDIR/${STRAIN}/${STRAIN}.orfs.fasta -evalue 0.00001 -outfmt 6 -max_hsps 1 -max_target_seqs 1 -out $OUT/${STRAIN}.CAP.out


done
