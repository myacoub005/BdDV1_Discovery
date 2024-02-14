#!/bin/bash
#SBATCH -p short -N 1 -n 16 --mem 16 --out logs/minimap.log

module load minimap2

# we need to make a paf file to plot linkages between the orthologs of scaffold_10 in other strains #
minimap2 -X -N 50 -p 0.1 -c scaffolds/Integration_scaffolds.fasta scaffolds/Integration_scaffolds.fasta > paf/Integration.paf
