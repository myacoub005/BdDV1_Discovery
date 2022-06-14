library(dplyr)
library(rvcheck)
library(ggplot2)
library(ggtree)
library(treeio)
library(ggstance)
library(ggtreeExtra)
library(RColorBrewer)
library(nationalparkcolors)

geneCopies <- read.table("virus_average_cov.txt", header=TRUE, sep="\t", row.names = NULL)
geneCopies

tip_metadata <- read.table("Lineages.tab", sep="\t", header=TRUE,check.names=FALSE, stringsAsFactor=F)

tip_metadata

tree <- read.tree("Bd_straintree_midpointroot")
tipnames <- tree$tip.label
tipnames
to_drop <- setdiff(tree$tip.label, geneCopies$Strain)
to_drop
straintree <- drop.tip(tree, to_drop)
tipnames <- straintree$tip.label

tipnames

p0 <- ggtree(tree, layout="circular") +
  geom_tiplab(size=0, color="black")

p0

p1 <- ggtree(straintree, layout="rectangular") +
  geom_tiplab(size=3, color="black")

p1

p <- p1 %<+% tip_metadata + geom_tippoint(aes(color=Lineage, shape=Continent), size=3) +
  scale_color_manual(values = park_palette("Everglades"))
plot(p)

p <- p1 %<+% tip_metadata + geom_tippoint(aes(color=Lineage), size=3) +
  scale_color_manual(values = park_palette("Everglades"))
plot(p)
#tip_metadata <- read.table("Lineages.tab", sep="\t", header=TRUE,check.names=FALSE, stringsAsFactor=F)

#tip_metadata

#p <- p1 %<+% tip_metadata + geom_tippoint(aes(color=Lineage), size=2)
#plot(p)

difftable <- setdiff(geneCopies$Strain, straintree$tip.label)
geneCopiesFilter <- filter(geneCopies,geneCopies$Strain %in% straintree$tip.label)

geneCopiesFilter
dd = data.frame(id=straintree$tip.label, value=(geneCopiesFilter$CAP))
dd

geneCounts = data.frame(geneCopiesFilter)

geneCounts

# Define the number of colors you want
#nb.cols <- 21
#mycolors <- colorRampPalette(brewer.pal(8, "Set3"))(nb.cols)
# Create a ggplot with 21 colors
# Use scale_fill_manual

ptbl <- facet_plot(p, panel = 'Vir_Ratio', data = geneCopiesFilter, geom = geom_barh, mapping = aes(x=VirRatio, group = label, fill=Lineage),
                   stat = "identity") + scale_fill_manual(values = park_palette("Everglades")) + theme_tree2(legend.position="bottom") 

ptbl
