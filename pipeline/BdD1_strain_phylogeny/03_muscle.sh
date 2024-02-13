#!/bin/bash
#SBATCH -p short -N 1 -n 16 --mem 16gb --out logs/muscle.log

module load muscle

#Combine the Cap nt sequences into one file and align

muscle -align BdDV1_CAP.fna -output BdDV1_CAP.fna.fasaln

module load iqtree

#iqtree2 -s reps.aa.fasaln
iqtree2 -s BdDV1_CAP.fna.fasaln -B 1000
