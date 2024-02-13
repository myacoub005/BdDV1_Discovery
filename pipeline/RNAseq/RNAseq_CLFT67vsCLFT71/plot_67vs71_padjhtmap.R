library(DESeq2)
library(tximport)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(magrittr)
library(Biobase)
library(pheatmap)
library(RColorBrewer)
library(ggplotify)
library(patchwork)


samples <- read.csv("strains_reduced.tsv",sep=",",header=TRUE)
samples <- subset(samples, select = -c(fq1, fq2) )
exprnames <- do.call(paste,c(samples[c("sample","unit")],sep="."))
files <- file.path("results",exprnames,"abundance.h5")
txi.kallisto <- tximport(files, type = "kallisto", txOut = TRUE)
head(txi.kallisto$counts)
colnames(txi.kallisto$counts) <- exprnames
colnames(txi.kallisto$abundance) <- exprnames
write.csv(txi.kallisto$abundance,"reports/kallisto_67_71.TPM.csv")
write.csv(txi.kallisto$counts,"reports/kallisto_67_71.counts.csv")

samples$Strain <- c("CLFT067", "CLFT067", "CLFT067", "CLFT071", "CLFT071", "CLFT071")
Strain = samples$Strain

# DEseq2 analyses


samples$BdDV1 <- c("Infected", "Infected", "Infected", "Uninfected", "Uninfected", "Uninfected")
BdDV1 = samples$BdDV1

sampleTable <- data.frame(BdDV1, Strain)
rownames(sampleTable) = exprnames

dds <- DESeqDataSetFromTximport(txi.kallisto,sampleTable, ~ BdDV1)

#nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 1, ]
#nrow(dds)

dds <- estimateSizeFactors(dds)
vsd <- vst(dds, blind=FALSE)
#rld <- rlog(dds, blind=FALSE)
head(assay(vsd), 3)

df <- bind_rows(
  as_tibble(log2(counts(dds, normalized=TRUE)[, 1:2]+1)) %>%
         mutate(transformation = "log2(x + 1)"),
 # as_tibble(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"),
  as_tibble(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"))

colnames(df)[1:2] <- c("x", "y")

#pdf("plots/RNASeq67_71_padj.pdf")
ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid( . ~ transformation)

#use DEGs 
res.sig.degs <- read.csv("reports/deseq_kallisto/Result_67_71_DEGs.csv")
select_mean <- order(res.sig.degs$padj)[1:50]
res.sig.degs.select.mean <- rownames(res.sig.degs[select_mean,])



#heatmap_data <- assay(vsd)[select,]
heatmap_data <- assay(vsd)[res.sig.degs.select.mean,]

colnames(heatmap_data) <- c("CLFT067.rep1", "CLFT067.rep2", "CLFT067.rep3", "CLFT071.rep1", "CLFT071.rep2", "CLFT071.rep3")
col_order <-c("CLFT067.rep1", "CLFT067.rep2", "CLFT067.rep3", "CLFT071.rep1", "CLFT071.rep2", "CLFT071.rep3")
heatmap_data <- heatmap_data[, col_order]

Strain <- c("CLFT067", "CLFT067", "CLFT067", "CLFT071", "CLFT071", "CLFT071")
BdDV1 <- c("Infected", "Infected", "Infected", "Uninfected", "Uninfected", "Uninfected")

sampd <- data.frame(BdDV1=BdDV1, Strain=Strain)
sampd#Strain <- factor(sampd$Strain, levels = c("CLFT067", "CLFT071"))

rownames(sampd) <- colnames(heatmap_data)

write.csv(heatmap_data, "reports/heatmap_50.csv")

ph <- as.ggplot(pheatmap(heatmap_data, cluster_rows=FALSE, show_rownames=TRUE,
         fontsize_row = 7,fontsize_col = 7,
         cluster_cols=FALSE, annotation_col=sampd,main=""))

ph

ggsave(filename = 'plots/heatmap_67_71_deg_top50_padj.pdf', plot = last_plot(), device = 'pdf', width = 7, height = 10, dpi = 300)
