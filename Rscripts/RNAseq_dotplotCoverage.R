#Reset R's Brain
rm(list=ls())
library(tidyverse)
library(purrr)
library(readr)
library(fs)
library(dplyr)
library(stringr)
library(ggplot2)
library(RColorBrewer)
library(scales)
library(viridis)
library(gggenomes)
library(patchwork)

dat <- read.table("BdDV1_44_coverage.txt", sep="\t", header=TRUE,check.names=FALSE, stringsAsFactor=F)
dat

dat$end <- dat$Start+1 #needs to have an "end" so +1
dat$file_id <- "coverage" #needs a "file id" or file name
dat <- dat %>%
  mutate(score = log(Coverage-4)/log(10))

dat$score
dat$Coverage
### log transform the average TPMs across sample 1:
### mutate(score = log(Coverage+1)/log(10))
Avg = log(96.0891)/log(10)
Avg
BT = log(1500)/log(10)
BT
#dat <- read_tsv("coverage.tab",
col_names = c("Scaffold","Start","Coverage") +
  col_types=c("ciid")
dat

plt1 <- ggplot(dat,aes(x = Blanked, y = score)) +
  geom_point() + scale_color_viridis(discrete=TRUE) + theme_bw() +
  labs(x="Position",y="log(Coverage+1)",title = "BdDV1 Coverage")
plt1

plt2 <- ggplot(dat,aes(x = Blanked, y = Coverage)) +
  geom_point() + scale_color_viridis(discrete=TRUE) + theme_bw() +
  labs(x="Position",y="Coverage",title = "BdDV1 Coverage")
plt2

library(cowplot)
plot_grid(plt1, plt2, labels = "AUTO")
# a minimal gene track
s0 <- tibble(
  seq_id = "BdDV1",
  length = 4405
)
g0 <- tibble(
  seq_id = c("BdDV1", "BdDV1", "BdDV1", "BdDV1"),
  Gene = c("REP", "CAP", "ORF3", "ORF4"),
  start = c(497, 1448, 2669, 3655),
  end = c(1411, 2316, 3375, 4278),
  strand = c("+","-","+","+"),
  introns = list(c(0,0), c(731,839), c(0,0), c(0,0)))
g0

p <- gggenomes(genes=g0, seqs=s0)
p

p1 <- p + 
  geom_seq(position=position_nudge(y=-1.5)) + # draw contig/chromosome lines
  geom_seq_label(nudge_y = -1.5) +   # label each sequence
  geom_gene(position=position_nudge(y=-1.5), aes(fill=seq_id)) + scale_fill_manual(values = c("#F8766D")) + theme(legend.position="left")
p1

px <- ggplot() + p1 +
  geom_line(data = dat, mapping = aes(x =dat$Blanked, y = score),position=position_nudge(y=0))
px
  
RNACov <- ggplot() + p1 +
  geom_segment(data = dat, mapping = aes(x =dat$Blanked,xend=(dat$Blanked +1), y = 0, yend=score),position=position_nudge(y=0)) +
  theme_bw() + labs(x="Position", y="log RNAseq Coverage",title = "BdDV1 RNAseq Coverage") + 
  geom_hline(yintercept = 1.98717, linetype=2, color = "#00BA38") + annotate("text", x = 2600, y=2.1, label = "Average Gene Expression", color = "#00BA38")  +
  geom_hline(yintercept = 3.176, linetype=2, color = "#0000FF") + annotate("text", x = 2600, y=3.3, label = "GADPH Expression", color = "#0000FF") +
  theme(legend.position="none", plot.background = element_blank())

RNACov

plot_grid(RNACov, labels = "AUTO")


