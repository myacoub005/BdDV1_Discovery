# install.packages("devtools")
#devtools::install_github("thackl/thacklr")
#devtools::install_github("thackl/gggenomes")

library(tidyverse)
library(gggenomes)
library(patchwork)
library(ggplot2)
library(gggenes)
library(ggfittext)
library(nationalparkcolors)

setwd("./results_for_plots/Integration")

#load in files, I subset to one scaffold before hand 
seq <- read_seqs("contig_138.fa") #fasta file
gff <- read_feats("contig_138.gff3") #gff

seq
gff

cov <- as.tibble(read.delim("VirExpr.txt", header = TRUE)) 
cov
cov$end <- cov$start+1 #needs to have an "end" so +1
cov$file_id <- "coverage" #needs a "file id" or file name
cov <- cov %>%
  mutate(score = log(tpm+1)) #mutating coverage data to be log normalized, added a pseudo count of +1 - log(1) = 0


#make gggenomes object and plot coverage info using geom_wiggle
p <- gggenomes(genes=g0, seqs=s0) +
  geom_wiggle(aes(z=score), fill="lavenderblush3", offset=-.3, height=.5, bounds=c(0,0,max(cov$score))) + 
  geom_seq() 
p
p <- gggenomes(seqs=seq, genes=gff)

p2 <- gggenomes(seq, gff) +
  geom_seq() + geom_bin_label() 

p2


# a minimal seq track
s0 <- tibble(
  seq_id = "BdDV1",
  length = 7000
)
s0
# a minimal gene track
g0 <- tibble(
  seq_id = c("BdDV1", "BdDV1", "BdDV1", "BdDV1"),
  gene_id = c("ORF1", "ORF2", "ORF3", "ORF4"),
  start = c(1, 904, 1963, 2894),
  end = c(624, 1527, 2831, 3742)
)


g0

p <- gggenomes(genes=g0, seqs=s0)
p + 
  geom_seq() +         # draw contig/chromosome lines
  geom_seq_label() +   # label each sequence 
  geom_gene() 

p1 <-p + gggenomes(seqs=s0, genes=g0) +
  geom_wiggle(aes(z=score), fill="lavenderblush3", offset=-.3, height=.5, bounds=c(0,0,max(cov$tpm))) + 
  geom_seq()
p1


q <- p + geom_feat(size=5,  position="identity", aes(color=g0$gene_id))
 
## Let's try an alignment with gggenomes 

#emale_seqs <- read.fasta("Virus_scaffold.fasta")
#emale_seqs

#emale_genes <- read_gff("Virus.gff")
#emale_genes
#g0
####### START AGAIN #####
seqID = c("scaffold_10","tig00000311")                  

s0 <- tibble(
  seq_id = c("CLFT044", "CLFT071"),
  length = c(972934,956479) 
)
s0

g0 <- read.table("int_features.txt", header = TRUE)
g0

tirs <- read.table("Virus_tirs.txt", header = TRUE)
tirs

p2 <- gggenomes(seqs = s0, genes = g0)
p2
p2 + 
  geom_seq() +         # draw contig/chromosome lines
  geom_seq_label() + 
  geom_gene(position = "strand") + geom_bin_label()

seqID = c("scaffold_10","tig00000311")


paf <- read.table("virus_Int.paf.txt", header = TRUE)
paf
p3 <- gggenomes(seqs = s0, genes = g0, links=paf, feats = tirs)
p3 +
geom_seq() + # draw contig/chromosome lines
geom_seq_label() +
geom_feat(size=5) +
geom_gene(position = "strand") +
geom_link()

p4 <- flip_by_links(p3, link_track = 1, min_coverage = 0.1)
global <- p4 +
geom_seq() + # draw contig/chromosome lines
geom_feat(size=15, color = "red") +
geom_gene(position = "strand") +
geom_link() + geom_gene(aes(fill=Organism, color = Organism),position = "strand") + coord_cartesian(ylim = c(0,5)) +
theme(legend.position="top", axis.text.y = element_blank()) + theme(axis.title.y=element_blank())
global


p5 <- p4 %>% flip(2)
p4 + 
  geom_seq() + 
  geom_feat(size=5) +# draw contig/chromosome lines
  geom_gene() + geom_link() + geom_gene(aes(fill=cluster_id),position = "strand")



example_genes
example_subgenes

ggplot(example_genes, aes(xmin = start, xmax = end, y = molecule, forward = orientation)) +
  facet_wrap(~ molecule, scales = "free", ncol = 1) +
  geom_gene_arrow(fill = "white") +
  geom_subgene_arrow(data = example_subgenes,
                     aes(xmin = start, xmax = end, y = molecule, fill = gene,
                         xsubmin = from, xsubmax = to), color="black", alpha=.7) +
  theme_genes()
genes <- read.table("Virus_genes.txt", header = TRUE)
subgenes <- read.table("Virus_subgenes.txt", header = TRUE)
genes
subgenes

ggplot(genes, aes(xmin = start, xmax = end, y = molecule, forward = orientation)) +
  facet_wrap(~ molecule, scales = "free", ncol = 1) +
  geom_gene_arrow(fill = "white") +
  geom_subgene_arrow(data = subgenes,
                     aes(xmin = start, xmax = end, y = molecule, fill = organism,
                         xsubmin = from, xsubmax = to), color="black", alpha=.7) +
  theme_genes()

detach("package:gggenomes", unload=TRUE)

genes <- read.table("Virus_int_closer.txt", header = TRUE)

genes

dummies <- make_alignment_dummies(
  genes,
  aes(xmin = start, xmax = end, y = molecule, id = gene),
  on = "CLFT044_004041"
)

dummies

ggplot(genes, aes(xmin = start, xmax = end, y = molecule, forward = orientation)) +
  facet_wrap(~ molecule, scales = "free", ncol = 1) +
  geom_gene_arrow(aes(fill=organism)) +
  theme_genes()

genes <- read.table("Virus_int_even_closer.txt", header = TRUE, fill = TRUE, sep = '\t')

genes

genes$gene2 <- gsub('(?=(?:.{4})+$)', "\n", genes$gene, perl = TRUE)

zoom <- ggplot(genes, aes(xmin = start, xmax = end, y = Scaffold_10, forward = orientation, label = locus)) +
  geom_gene_arrow(aes(fill=organism, color =organism)) + scale_fill_manual(values = c("#F8766D", "#00BA38", "#F8766D")) +
  geom_text(data=genes %>% mutate(start = (start + end)/2), aes(x=start, label = locus, angle =90), nudge_y = 9.5 ) +
  coord_cartesian(ylim = c(0,20)) +
  theme(legend.position="none", axis.text.y = element_blank()) + labs(x="Location (bp)")

zoom

library(cowplot)

plot_grid(global, zoom, labels = "AUTO", ncol = 1)



### Plot integration locus Coverage for a set of strains ####

# First we need to make a subsetted scaffold and gene set to compare against 
seqID = c("scaffold_10")

s0 <- tibble(
  seq_id = c("scaffold_10"),
  length = c(150086) 
)
s0

g0 <- read.table("Virus_int_closer.txt", header = TRUE)
g0

p <- gggenomes(seqs = s0, genes = g0)
p
p + 
  geom_seq() +         # draw contig/chromosome lines
  geom_seq_label() + 
  geom_gene(position = "strand")

closer <- p +
  geom_seq() + # draw contig/chromosome lines
  geom_gene(position = "strand") +
  geom_gene(aes(fill=Organism),position = "strand") +
  theme(legend.position="bottom", axis.text.y = element_blank()) + theme(axis.title.y=element_blank())

closer

cov <- as.tibble(read.delim("scaffold_10_coverage.txt", header = TRUE))
cov

p2 <- closer + geom_segment(data=cov, aes(x=start, xend=end, y=Coverage.Average_NCL, yend=Coverage.Average_NCL, color=Lineage, linetype = BdDV1, size=5), position=position_nudge(y=1.5))
p2  

p2 <- closer + geom_line(data=cov, aes(x=start, xend=end, y=Coverage.Average_NCL, yend=Coverage.Average_NCL, color=Lineage, linetype = BdDV1), position=position_nudge(y=1.5)) + scale_color_manual(values = park_palette("Everglades"))
p2  

### Let's take a closer look ####

# First we need to make a subsetted scaffold and gene set to compare against 
seqID = c("scaffold_10")

s0 <- tibble(
  seq_id = c("scaffold_10"),
  length = c(36500) 
)
s0

g0 <- read.table("Virus_coverage_close.txt", header = TRUE)
g0

p <- gggenomes(seqs = s0, genes = g0)
p
p + 
  geom_seq() +         # draw contig/chromosome lines
  geom_seq_label() + 
  geom_gene(position = "strand")

closer <- p +
  geom_seq() + # draw contig/chromosome lines
  geom_gene(position = "strand") +
  geom_gene(aes(fill=Organism),position = "strand") +
  theme(legend.position="bottom", axis.text.y = element_blank()) + theme(axis.title.y=element_blank())

closer

cov <- as.tibble(read.delim("scaffold_10_coverage_closer.txt", header = TRUE))
cov

p2 <- closer + geom_segment(data=cov, aes(x=start, xend=end, y=Coverage.Average_NCL, yend=Coverage.Average_NCL, color=Lineage, linetype = BdDV1), position=position_nudge(y=1.25))
p2 

p2 <- closer + geom_line(data=cov, aes(x=start, y=Coverage.Average_NCL, color=Lineage, group=Strain, linetype = BdDV1), position=position_nudge(y=1)) + scale_color_manual(values = park_palette("Everglades"))
p2


# First we need to make a subsetted scaffold and gene set to compare against 
seqID = c("scaffold_10")

s0 <- tibble(
  seq_id = c("scaffold_10"),
  length = c(972934) 
)
s0

g0 <- read.table("scaffold_10_loci.txt", header = TRUE)
g0

p <- gggenomes(seqs = s0, genes = g0)
p
p + 
  geom_seq() +         # draw contig/chromosome lines
  geom_seq_label() + 
  geom_gene(position = "strand")

closer <- p +
  geom_seq() + # draw contig/chromosome lines
  geom_gene(position = "strand") +
  geom_gene(aes(fill=Organism),position = "strand") +
  theme(legend.position="bottom", axis.text.y = element_blank()) + theme(axis.title.y=element_blank())

closer

cov <- as.tibble(read.delim("scaffold_10_total_coverage.txt", header = TRUE))
cov

p2 <- closer + geom_segment(data=cov, aes(x=start, xend=end, y=Coverage.AverageNCL, yend=Coverage.AverageNCL, color=Lineage, linetype = BddV1), position=position_nudge(y=1.25))
p2 

p2 <- closer + geom_line(data=cov, aes(x=start, y=Coverage.AverageNCL, color=Lineage, group=Strain, linetype = BddV1), position=position_nudge(y=1)) + scale_color_manual(values = park_palette("Everglades")) + scale_y_continuous(limits = c(0,2.5)) +
  coord_cartesian(ylim = c(0,5))
p2

values <- loess(cov$Coverage.AverageNCL ~ cov$start)
values

###### And let's start yet again but this time with JEL423 FFS ###

####### START AGAIN #####
              

s0 <- tibble(
  seq_id = c("JEL423", "CLFT044", "CLFT071"),
  length = c(972934,956479,979369) 
)
s0

g0 <- read.table("int_features.JEL423.txt", header = TRUE)
g0

tirs <- read.table("Virus_tirs.txt", header = TRUE)
tirs

p2 <- gggenomes(seqs = s0, genes = g0)
p2
p2 + 
  geom_seq() +         # draw contig/chromosome lines
  geom_seq_label() + 
  geom_gene(position = "strand") + geom_bin_label()

seqID = c("scaffold_10","tig00000311","DS022310.1")


paf <- read.table("virus_Int.JEL423.paf.txt", header = TRUE)
paf
p3 <- gggenomes(seqs = s0, genes = g0, links=paf, feats = tirs)
p3 +
  geom_seq() + # draw contig/chromosome lines
  geom_seq_label() +
  geom_feat(size=5) +
  geom_gene(position = "strand") +
  geom_link()

p3

p4 <- flip_by_links(p3, link_track = 1, min_coverage = 0.1)

p5 <- p3 %>% flip(3)

global <- p5 +
  geom_seq() + # draw contig/chromosome lines
  geom_feat(size=15, color = "red") +
  geom_gene(position = "strand") +
  geom_link() + geom_gene(aes(fill=Organism, color = Organism),position = "strand") +
  scale_fill_manual(values = c("#F8766D", "#00BA38", "#3399FF", "#9966FF")) +
  scale_color_manual(values = c("#F8766D", "#00BA38", "#3399FF", "#9966FF")) +
  theme(legend.position="top", axis.text.y = element_blank()) + theme(axis.title.y=element_blank())
global



example_genes
example_subgenes

ggplot(example_genes, aes(xmin = start, xmax = end, y = molecule, forward = orientation)) +
  facet_wrap(~ molecule, scales = "free", ncol = 1) +
  geom_gene_arrow(fill = "white") +
  geom_subgene_arrow(data = example_subgenes,
                     aes(xmin = start, xmax = end, y = molecule, fill = gene,
                         xsubmin = from, xsubmax = to), color="black", alpha=.7) +
  theme_genes()
genes <- read.table("Virus_genes.txt", header = TRUE)
subgenes <- read.table("Virus_subgenes.txt", header = TRUE)
genes
subgenes

ggplot(genes, aes(xmin = start, xmax = end, y = organism, forward = orientation)) +
  facet_wrap(~ organism, scales = "free", ncol = 1) +
  geom_gene_arrow(fill = "white") +
  geom_subgene_arrow(data = subgenes,
                     aes(xmin = start, xmax = end, y = organism, fill = organism,
                         xsubmin = from, xsubmax = to), color="black", alpha=.7) +
  theme_genes()
