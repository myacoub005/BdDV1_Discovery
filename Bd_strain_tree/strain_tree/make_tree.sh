#!/usr/bin/bash -l
#SBATCH -p intel -N 1 -n 6 --mem 64gb 
module load iqtree

iqtree2 -s Bd_2.All.SNP.mfa -m GTR+ASC -st DNA -B 1000 --alrt 1000 -nt AUTO

