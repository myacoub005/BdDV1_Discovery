#!/bin/bash
#SBATCH -p short -N 1 -n 16 --mem 16gb --out logs/muscle.log

module load muscle

muscle -align REP.CRESS.fasta -output REP.CRESS.fasaln

