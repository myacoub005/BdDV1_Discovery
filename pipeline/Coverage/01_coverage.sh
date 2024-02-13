#!/usr/bin/bash

#SBATCH --nodes 1 --ntasks 24 --mem 24G -p short -J readcov --out logs/bwa.%a.log --time 2:00:00

#Assemblies are located here /rhome/myaco005/shared/projects/Chytrid/Bd_popgen/Asm_Unmapped_BdVir

#ln -s /rhome/myaco005/shared/projects/Chytrid/Bd_popgen/Asm_Unmapped_BdVir/input

module unload miniconda2
module load miniconda3
module load samtools
module load bwa
module load bam2fastq
module load mosdepth

CPU=$SLURM_CPUS_ON_NODE
N=${SLURM_ARRAY_TASK_ID}
TEMP=/scratch
if [ ! $N ]; then
    N=$1
    if [ ! $N ]; then
        echo "Need an array id or cmdline val for the job"
        exit
    fi
fi

INDIR=input
OUT=ref_aln
COV=coverage
SAMPLEFILE=samples_PE.dat
REF=genome/GCA_000149865.1_BD_JEL423_genomic_BdDV1.fa
REFBED=genome/chroms.bed
if [ ! -f $REF.bwt ]; then
    bwa index $REF
fi

mkdir -p $COV $OUT
sed -n ${N}p $SAMPLEFILE | while read BASE STRAIN
do
    BAM=$OUT/$STRAIN.unmapped.bam
    ALLBAM=$OUT/$STRAIN.all.cram
    if [ ! -f $COV/${STRAIN}.regions.bed.gz ]; then
	if [ ! -f $ALLBAM ]; then
	    bwa mem -t $CPU $REF $INDIR/${BASE}_[12].fastq.gz | samtools sort --threads $CPU --reference $REF -Ocram -o $ALLBAM -
	    samtools index $ALLBAM
	    samtools view -h --threads $CPU -u -f4 -Obam -o $BAM $ALLBAM
	fi
	mosdepth -f $REF -x -n --by $REFBED -t $CPU $COV/$STRAIN $ALLBAM
fi
    if [ ! -f $OUT/${STRAIN}_1.fastq.gz ]; then    
	bam2fastq --no-aligned -o $OUT/${STRAIN}#.fastq $BAM
	find $OUT -name "${STRAIN}_*" -size 0 | xargs rm
	pigz $OUT/${STRAIN}_[12].fastq
    fi
done
