#!/usr/bin/bash
#SBATCH -p short -N 1 -n 8 --mem 8gb --out logs/kallisto.%a.log

CPU=$SLURM_CPUS_ON_NODE

if [ ! $CPUS ]; then
    CPU=1
fi

N=${SLURM_ARRAY_TASK_ID}

if [ ! $N ]; then
    N=$1
    if [ ! $N ]; then
        echo "need to provide a number by --array or cmdline"
        exit
    fi
fi

module load kallisto
DB=genome/JEL423_CLFT044.idx
OUTDIR=$(realpath results)
INDIR=$(realpath fastq)

IFS=,
tail -n +2 strains.updated.tsv | sed -n ${N}p | while read SAMPLE REP FQ1 FQ2
do
	OUTNAME=$SAMPLE.$REP
        kallisto quant -i $DB -t $CPU -o $OUTDIR/$OUTNAME --bias $FQ1 $FQ2
done
