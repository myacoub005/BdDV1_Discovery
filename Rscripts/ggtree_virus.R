library(dplyr)
library(rvcheck)
library(ggplot2)
library(ggdendro)
library(phangorn)
library(ggtree)
library(treeio)
library(ggstance)
library(ggtreeExtra)
library(RColorBrewer)
library(nationalparkcolors)
library(ggarange)
library(gridGraphics)
library(patchwork)
library("ggplotify")
library("grid")

setwd("results_for_plots/Viral_Presence")

geneCopies <- read.table("virus_average_cov.txt", header=TRUE, sep="\t", row.names = NULL)
geneCopies

tip_metadata <- read.table("Lineages.tab", sep="\t", header=TRUE,check.names=FALSE, stringsAsFactor=F)

tip_metadata

#tree <- read.tree("Bd_straintree_midpointroot")
tree <- read.tree("Bd_straintree_Jason")
tipnames <- tree$tip.label
tipnames
to_drop <- setdiff(tree$tip.label, geneCopies$Strain)
to_drop
straintree <- drop.tip(tree, to_drop)
tipnames <- straintree$tip.label

tipnames

p0 <- ggtree(tree, layout="rectangular") +
  geom_tiplab(size=3, color="black")

p0

p1 <- ggtree(straintree, layout="rectangular") +
  geom_tiplab(size=5, color="black", offset=0.01)

p1

p <- p1 %<+% tip_metadata + geom_tippoint(aes(color=Lineage, shape=Continent), size=2.5) +
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

ptbl <- facet_plot(p, panel = 'Virus Ratio', data = geneCopiesFilter, geom = geom_barh, mapping = aes(x=VirRatio, group = label, fill=Lineage),
                   stat = "identity") + scale_fill_manual(values = park_palette("Everglades")) + theme_tree2(legend.position="none", axis.text.x =element_text(size=12, face="bold")) + xlim_expand(c(0,0.4), panel = "Tree")

ptbl

geneCopiesFilter


### Plot the Overall simplified Barplots #####

Presence <- read.table("Simplified_BdDV1_coverage.txt", header = TRUE)
Presence

bar <- ggplot(data = Presence, aes(x=Lineage, y=Percent_Virus., fill = Lineage)) + geom_bar(stat="identity") +
  scale_fill_manual(values = park_palette("Everglades")) + geom_text(aes(label = N), vjust=-0.3, size=10) + theme_minimal() + ylab("Percent Virus +") + theme(legend.position = "left", axis.title.x = element_text(size = 0),axis.title.y = element_text(size = 16),axis.text.y =element_text(size=12, face="bold"),axis.text=element_text(size=12, face="bold"), legend.key.size =unit(1.5, 'cm'), legend.title = element_text(size=14), legend.text = element_text(size=12))

bar
###### Can we do a co-phylogeny sorta thing? ######

library(ape)
library(phytools)
library(dendextend)
library(viridis)
library(dplyr)
library(phylogram)

#BdVirus.CAP.tree worked better
#tree2 <- read.tree("BdVirus.CAP.tree")
tree2 <- read.tree("BdDV1_CAP")
tree1 <- tree
tipnames <- tree2$tip.label
tipnames
#drop tips not in tree2
to_drop <- setdiff(tree1$tip.label, tree2$tip.label)
to_drop
tree1 <- drop.tip(tree1, to_drop)

#drop tips not in tree1
to_drop <- setdiff(tree2$tip.label, tree1$tip.label)
to_drop
tree2 <- drop.tip(tree2, to_drop)


tree1 <- compute.brlen(tree1)
tree2 <- compute.brlen(tree2)
tree1<- as.dendrogram(tree1)
tree1
tree2<- as.dendrogram(tree2)
tree2


dndlist <- dendextend::dendlist(tree1, dend15_2)
dendextend::tanglegram(dndlist, fast = TRUE, margin_inner = 1, lab.cex = 2, lwd = 
                         0.5, edge.lwd = 0.75, type = "r")
dev.copy(pdf, 'Discrete001.pdf', width = 10, height = 11)
dev.off()


dend15_2 <- tree2 %>% 
  set("labels_col", c(2)) %>%  # change color 
  set("labels_cex", c(2))
dend15_2 %>% plot(main = "After")

dend12 <- dendlist(tree1, tree2)
dev.new(width=5, height=4)
dend12 %>% tanglegram(common_subtrees_color_branches = TRUE,
                      lab.cex = .5, margin_inner = 2.3)
dev.copy(pdf, 'BdVirus_tangelgram.pdf', width = 10, height = 11)

dend12 %>% tanglegram(common_subtrees_color_branches = TRUE, axes = FALSE, 
                      lab.cex = .5, margin_inner = 2.3)

labels <- tree2 %>% set("labels_to_char") %>% labels 
labels <- as.data.frame(labels)
labels
labels2 <- merge(labels, tip_metadata, by.x="labels", by.y="Strain", sort=F)
labels2
cols <- as.character(labels2$Color)
cols

colors <- c("Green", "Green", "Green", "Green", "Green", "Red", "Red","Red","Red","Blue","Blue","Blue","Blue","Blue")

#CAP gene color code
#colors <- c("Green", "Green","Green","Green","Green","Green","Green","Green","Green","Green","Green","Green","Green","Green","Green","Green","Green","Green","Green","Red","Red","Red","Red","Red","Red","Red","Red","Blue","Blue","Blue","Blue","Blue","Blue","Blue","Blue","Blue")

colors <- c("Sienna", "Sienna","Sienna","Sienna","Sienna","Sienna","Sienna","Sienna","Sienna","Sienna","Sienna","Sienna","Sienna","Sienna","Sienna","Sienna","Sienna","Sienna","darkgoldenrod1","azure4","darkcyan","darkcyan","darkcyan","darkcyan","darkcyan","darkcyan","darkcyan","cadetblue1","cadetblue1","cadetblue1","cadetblue1","cadetblue1","cadetblue1","cadetblue1","cadetblue1","cadetblue1")


tanglegram(dend12, color_lines = colors)

# Let's find the optimal entanglement
set.seed(3958)
x <- dend12 %>% untangle(method = "random", R = 10) 
x %>% plot(main = paste("entanglement =", round(entanglement(x), 2)))

dend12 %>% 
  set("leaves_pch", 21) %>% 
  set("leaves_bg", "gold") %>%    
  set("leaves_cex", 2) %>% 
  set("leaves_col", "darkred") %>% 
  plot(main = "(6) Show (larger+colored+filled)\n leaves")

x <- dend12 %>% untangle(method = "step2side") 
x %>% plot(main = paste("entanglement =", round(entanglement(x), 2)))
x %>% set("labels_cex", 30)
#colors <- c("Red", "Red","Red","Red","Red","Red","Red","Red","Red", "Red","Red","Red","Red","Red","Red","Red","Red", "Red","Green","Green","Green","Green","Green","Green","Green","Green","Blue","Blue","Blue","Blue","Blue", "Blue", "Blue")

# Plot tanglegram of CAP vs Bd + color lines by lineage
tangle <- tanglegram(x, color_lines = colors, lab.cex=1.9, margin_inner = 12, main_left="Bd Phylogeny", main_right ="BdDV1 Phylogeny", match_order_by_labels=TRUE)
tangle

ggplot(tangle, theme = NULL) # horiz plot (and let's remove theme) in ggplot2

# Calculate Correlation

dist.dendlist(dend12)
cor.dendlist(dend12)

hcdata1 <- dendro_data(tree1, 5, type = "rectangle")
hcdata2 <- dendro_data(tree2, 5, type = "rectangle")






########### test

geneCopies <- read.table("virus_average_cov.txt", header=TRUE, sep="\t", row.names = NULL)
geneCopies

tip_metadata <- read.table("Lineages.tab", sep="\t", header=TRUE,check.names=FALSE, stringsAsFactor=F)

tip_metadata

#tree <- read.tree("Bd_straintree_midpointroot")
tree <- read.tree("Bd_pang_tree_bs_Jason.txt")
tree$node.label

tree1 <- ggtree(tree) + geom_text2(aes(subset=!isTip, label=node), hjust=-.3) + geom_tiplab()
tree1
tipnames <- tree$tip.label
tipnames
to_drop <- setdiff(tree$tip.label, geneCopies$Strain)
to_drop
straintree <- drop.tip(tree, to_drop)
tipnames <- straintree$tip.label

tipnames

p0 <- ggtree(tree, layout="rectangular") +
  geom_tiplab(size=3, color="black")

p0

p1 <- ggtree(straintree, layout="rectangular") +
  geom_tiplab(size=5, color="black", offset=0.01)

p1

p <- p1 %<+% tip_metadata + geom_tippoint(aes(color=Lineage, shape=Continent), size=2.5) +
  scale_color_manual(values = park_palette("Everglades"))
plot(p)

p <- p1 %<+% tip_metadata + geom_tippoint(aes(color=Lineage), size=3) +
  scale_color_manual(values = park_palette("Everglades"))
plot(p)


p1 <- ggtree(straintree, layout="rectangular") +
  geom_nodelab(label= straintree$node.label,size=2, node = "internal") +
  geom_tiplab(size=3.5, color="black", offset=0.01) 
p1


p <- p1 %<+% tip_metadata + geom_tippoint(aes(color=Lineage), size=3) +
  geom_nodelab(label= straintree$node.label,size=3, node = "internal") +
  scale_color_manual(values = park_palette("Everglades"))
plot(p)


ptbl <- facet_plot(p, panel = 'Virus Ratio', data = geneCopiesFilter, geom = geom_barh, mapping = aes(x=VirRatio, group = label, fill=Lineage),
                   stat = "identity") + scale_fill_manual(values = park_palette("Everglades")) + theme_tree2(legend.position="none", axis.text.x =element_text(size=12, face="bold")) + xlim_expand(c(0,0.4), panel = "Tree")

ptbl
