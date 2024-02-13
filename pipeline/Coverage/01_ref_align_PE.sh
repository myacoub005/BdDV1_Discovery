#!/usr/bin/bash
#SBATCH -p short -n 24 --mem 64gb --out logs/bwa_PE.%a.log -N 1

module load samtools
module load bwa

CPU=$SLURM_CPUS_ON_NODE
N=${SLURM_ARRAY_TASK_ID}

if [ ! $N ]; then
    N=$1
    if [ ! $N ]; then
        echo "Need an array id or cmdline val for the job"
        exit
    fi
fi

INDIR=input
OUT=aln_TF5a1
SAMPLEFILE=samples_PE.dat
TEMP=/scratch
#REF=genome/assembled_TF5a1.fa
REF=assembled_TF5a1.fa

if [ ! -f $REF.bwt ]; then
    bwa index $REF
fi

mkdir -p $OUT

sed -n ${N}p $SAMPLEFILE | while read LIBRARY STRAIN
do
    ALNFILE=$OUT/$STRAIN.cram
    echo "$LIBRARY -> $ALNFILE"
    if [ ! -s $ALNFILE ]; then
	    bwa mem -t $CPU $REF $INDIR/${LIBRARY}_[12].fastq.gz | samtools fixmate --threads $CPU -Ocram --reference $REF - $TEMP/${LIBRARY}.fix.cram
	    samtools sort --threads $CPU --reference $REF -Ocram -o $ALNFILE $TEMP/${LIBRARY}.fix.cram
samtools index $ALNFILE
#samtools view -h -b $ALNFILE $TEMP/${LIBRARY}.fix.cram "contig_138:3537-7941" > $TEMP/${LIBRARY}.bam
#samtools index $TEMP/${LIBRARY}.bam
samtools fastq $ALNFILE > Virus_fq/$STRAIN.fq
#rm $TEMP/${LIBRARY}.*
    else
	echo "Not running already have $ALNFILE"
    fi
done


