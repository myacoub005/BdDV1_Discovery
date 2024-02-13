#!/bin/bash
#SBATCH -p intel
#SBATCH --time 04:00:00
#SBATCH --ntasks 10
#SBATCH --nodes 1
#SBATCH --mem 48G
#SBATCH --out logs/iqtree.rep.log

module load iqtree

#iqtree2 -s reps.aa.fasaln
iqtree2 -s REP.CRESS.fasaln -B 1000
