#load libraries
library(DESeq2)
library(tximport)
library(dplyr)
library(ggplot2)
library(magrittr)
library(Biobase)
library(pheatmap)
library(RColorBrewer)

#read data in 

samples <- read.csv("strains_reduced.tsv",sep=",",header=TRUE)
samples <- subset(samples, select = -c(fq1, fq2) )
exprnames <- do.call(paste,c(samples[c("sample","unit")],sep="."))
files <- file.path("results",exprnames,"abundance.h5")
txi.kallisto <- tximport(files, type = "kallisto", txOut = TRUE)
head(txi.kallisto$counts)
colnames(txi.kallisto$counts) <- samples$Name
colnames(txi.kallisto$abundance) <- samples$Name

#going to group data here - otherwise would need to deal with contrasts and loops 
#just testing pipeline now so this is easier to start with
#plus big differences in heatmap from kallisto_DESeq_1var.R report show differences
#between these two groups
samples$BdDV1 <- c("Infected", "Infected", "Infected", "Uninfected", "Uninfected", "Uninfected")

samples$Strain <- c("CLFT067", "CLFT067", "CLFT067", "CLFT071", "CLFT071", "CLFT071")

# DEseq2 set-up
treatment = samples$BdDV1

sampleTable <- data.frame(condition=treatment)
sampleTable$condition <- factor(sampleTable$condition)

rownames(sampleTable) = samples$Name

#import data into deseq looking at "condition" which is Infected vs. Uninfected
dds <- DESeqDataSetFromTximport(txi.kallisto,sampleTable, ~ condition)

dds <- dds[ rowSums(counts(dds)) > 1, ]


dds <- estimateSizeFactors(dds)

#run deseq diff abundance analysis
dds.test <- DESeq(dds)

#get results
#this is where we would need contrasts / loops
#if we had more than two groups to compare
res <- results(dds.test) 

#get only results that are sig & have log fold change >2
res.sig <- subset(res, (padj < 0.05 & abs(log2FoldChange) >2))

#write table with all differentially expressed genes (DEG) between Infected and Uninfected
write.table(res, "reports/deseq_kallisto/All_67_71_genes.csv",
            row.names=T,
            sep=",")
