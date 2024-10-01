#define packages to install
packages <- c('ape')

#install all packages that are not already installed
install.packages(setdiff(packages, rownames(installed.packages())))

library(tidyverse)
library(ape)

args <- commandArgs(trailingOnly = TRUE)
treefile <- args[1]
outgroup <- args[2]

#treefile = 'RAxML_bipartitions.alignment_pi_m1_nogap'
#treefile='alignment_pi_m1_nogap.cf.tree'
#outgroup = 'AotNan'

tree <- read.tree(treefile)
tree <- root.phylo(tree, outgroup)

pdf(file = paste0(treefile,'.pdf'))
plot(tree, show.node.label=TRUE)
dev.off()