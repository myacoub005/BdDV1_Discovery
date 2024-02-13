#!/bin/bash
#SBATCH -p short -N 1 -n 2 --out logs/00_initialize.log

DB=genome
mkdir -p $DB

#download the ref genome for analysis

#curl -L https://ftp.ncbi.nlm.nih.gov/genomes/genbank/fungi/Batrachochytrium_dendrobatidis/latest_assembly_versions/GCA_000149865.1_BD_JEL423/GCA_000149865.1_BD_JEL423_cds_from_genomic.fna.gz | pigz -dc > $DB/GCA_000149865.1_BD_JEL423.fasta

#Now add the viral genes to the ref genome
cat ~/shared/projects/Chytrid/BdVirus/assembled_TF5a1.cds >> $DB/GCA_000149865.1_BD_JEL423.fasta

module load kallisto

kallisto index -i $DB/GCA_000149865.1.idx genome/GCA_000149865.1_BD_JEL423.fasta -k 27

ln -s /bigdata/stajichlab/shared/projects/SeqData/Novogene_SeqData/BdDV1_RNASeq/usftp21.novogene.com/raw_data
<<<<<<< HEAD
=======

mkdir -p fastq

rsync -a raw_data/*/*.gz fastq/
>>>>>>> c7827e723ef8caca6bf8ea5d9ca8a892815aef18

mkdir -p fastq

rsync -a raw_data/*/*.gz fastq/
