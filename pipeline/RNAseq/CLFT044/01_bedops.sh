#!/bin/bash
#SBATCH -p short -N 1 -n 16 --mem 16gb --out logs/bedops.log

module load bedops

grep "scaffold_10" coverage.res.tsv > scaffold_10_coverage.tsv
awk -vFS=" " -vOFS=" " '{ print $1, $2, ($2 + 1), ".", $3 }' scaffold_10_coverage.tsv > test.bed
bedops --merge test.bed | bedops --chop 10000 - | bedmap --echo --mean --delim ' ' - test.bed > answer.bed
