####This script will be for visualizing GAG presence across fungi ####

rm(list=ls())
install.packages("ggrepel")
library(dplyr)
library(rvcheck)
library(ggplot2)
library(ggtree)
library(treeio)
library(phytools)
library(ggstance)
library(ggtreeExtra)
library(RColorBrewer)
library(ggnewscale)
library(ggrepel)

tip_metadata <- read.table("CRESS_metadata.tab", sep="\t", header=TRUE,check.names=FALSE, stringsAsFactor=F)
geneCopies <- read.table("CRESS_metadata.tab", header=TRUE, sep="\t", row.names = 1)
geneCopies
tip_metadata

tree <- read.tree("Cress.tree")
tipnames <- tree$tip.label
tipnames
to_drop <- setdiff(tree$tip.label, rownames(geneCopies))
to_drop
straintree <- drop.tip(tree, to_drop)
tipnames <- straintree$tip.label
tipnames

p1 <- ggtree(straintree, layout="circular", branch.length = 'none') +
  geom_text2(aes(subset=!isTip, label=node)) +
  geom_tiplab(size=2, color="black")

p1 

tree1 <- ggtree(straintree) + geom_text2(aes(subset=!isTip, label=node), hjust=-.3) + geom_tiplab()
tree1
tree2 <- phytools::reroot(straintree, 98)
plot(tree2)
nodelabels()
title("phytools::reroot")

p1 <- ggtree(tree2, layout="circular") +
  geom_tiplab(size=2, color="black")
p1

p2 <- ggtree(tree2, branch.length = 'none', layout="fan") + 
  geom_cladelab(node=138, label="Geminiviridae", angle=310, 
                fontsize=5, offset=4.5, vjust=-2, hjust=0.5) + 
  geom_cladelab(node=111, label="Genomoviridae", angle=380, 
                fontsize=5, offset=4.5, vjust=2.5, hjust=0.5) +
  geom_cladelab(node=82, label="Circoviridae", angle=60, 
                fontsize=5, offset=4.65, vjust=-1.5, hjust=0) 
p <- p2 %<+% tip_metadata + geom_tippoint(aes(color=Host), size=3) +
  scale_color_manual(values = c("skyblue4", "red3", "forestgreen", "goldenrod3", "plum4")) +
  geom_tiplab(size=2, color="black", offset = 0.5)

plot(p)


#####Let's try to add in the node support values ######
tree <- read.tree("CRESS.itol.tree.txt")
tree
##re-root the tree at geminiviridae clade
tree2 <- phytools::reroot(tree, 121)

tipnames <- tree2$tip.label
tipnames
tip_metadata <- read.table("CRESS_metadata.tab", sep="\t", header=TRUE,check.names=FALSE, stringsAsFactor=F)
tip_metadata
to_drop <- setdiff(tree2$tip.label, rownames(geneCopies))
to_drop
straintree <- drop.tip(tree2, to_drop)
tipnames <- straintree$tip.label
tipnames

p1 <- ggtree(straintree, layout="circular") +
  geom_nodelab(size=2, vjust = -0.5) +
  geom_tiplab(size=2, color="black", offset = 1) + 
  #geom_text2(aes(subset=!isTip, label=node)) +
  geom_cladelab(node=134, label="Geminiviridae", angle=-40, fontsize=5, offset.text=1, vjsut=-5, hjust=0.5, offset=15) +
  geom_cladelab(node=79, label="Circoviridae", angle=70, fontsize=5, offset.text=1, vjsut=-1.5, hjust=0.5, offset=15) +
  geom_cladelab(node=105, label="Genomoviridae", angle=380, fontsize=5, offset.text=1, vjsut=2.5, hjust=0.5, offset=15)

p1

p <- p1 %<+% tip_metadata + geom_tippoint(aes(color=Host), size=3) +
  scale_color_manual(values = c("skyblue4", "red3", "forestgreen", "goldenrod3", "plum4"))

plot(p)

p1 <- ggtree(straintree, layout="rectangular") +
  geom_nodelab(size=2) +
  geom_tiplab(size=3, color="black", offset=.1) + 
  #geom_text2(aes(subset=!isTip, label=node)) +
  geom_cladelab(node=134, label="Geminiviridae", fontsize=3, offset.text=0.5, vjsut=-5, hjust=0.5, offset=2) +
  geom_cladelab(node=79, label="Circoviridae", fontsize=3, offset.text=0.5, vjsut=-1.5, hjust=0.5, offset=1.2) +
  geom_cladelab(node=105, label="Genomoviridae", fontsize=3, offset.text=0.5, vjsut=2.5, hjust=0.5, offset=2)

p1

p <- p1 %<+% tip_metadata + geom_tippoint(aes(color=Host), size=3) +
  scale_color_manual(values = c("skyblue4", "red3", "forestgreen", "goldenrod3", "plum4")) + theme_tree2(legend.position="bottom")

plot(p)


#####Let's try to add in the node support values ######
tree <- read.tree("CRESS.itol.tree.txt")
tree
##re-root the tree at geminiviridae clade
tree2 <- phytools::reroot(tree, 121)

tipnames <- tree2$tip.label
tipnames
tip_metadata <- read.table("CRESS_metadata.tab", sep="\t", header=TRUE,check.names=FALSE, stringsAsFactor=F)
tip_metadata
to_drop <- setdiff(tree2$tip.label, rownames(geneCopies))
to_drop
straintree <- drop.tip(tree2, to_drop)
tipnames <- straintree$tip.label
tipnames

p1 <- ggtree(straintree, layout="circular") +
  geom_nodelab(size=2, vjust = -0.5) +
  geom_tiplab(size=2, color="black", offset = 1) + 
  #geom_text2(aes(subset=!isTip, label=node)) +
  geom_cladelab(node=134, label="Geminiviridae", angle=-40, fontsize=5, offset.text=1, vjsut=-5, hjust=0.5, offset=15) +
  geom_cladelab(node=79, label="Circoviridae", angle=70, fontsize=5, offset.text=1, vjsut=-1.5, hjust=0.5, offset=15) +
  geom_cladelab(node=105, label="Genomoviridae", angle=380, fontsize=5, offset.text=1, vjsut=2.5, hjust=0.5, offset=15)

p1

p <- p1 %<+% tip_metadata + geom_tippoint(aes(color=Host), size=3) +
  scale_color_manual(values = c("skyblue4", "red3", "forestgreen", "goldenrod3", "plum4"))

plot(p)

p1 <- ggtree(straintree, layout="rectangular") +
  geom_nodelab(size=2) +
  geom_tiplab(size=3, color="black", offset=.1) + 
  #geom_text2(aes(subset=!isTip, label=node)) +
  geom_cladelab(node=134, label="Geminiviridae", fontsize=3, offset.text=0.5, vjsut=-5, hjust=0.5, offset=2) +
  geom_cladelab(node=79, label="Circoviridae", fontsize=3, offset.text=0.5, vjsut=-1.5, hjust=0.5, offset=1.2) +
  geom_cladelab(node=105, label="Genomoviridae", fontsize=3, offset.text=0.5, vjsut=2.5, hjust=0.5, offset=2)

p1

p <- p1 %<+% tip_metadata + geom_tippoint(aes(color=Host), size=3) +
  scale_color_manual(values = c("skyblue4", "red3", "forestgreen", "goldenrod3", "plum4")) + theme_tree2(legend.position="bottom")

plot(p)
