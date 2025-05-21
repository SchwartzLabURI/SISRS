if (!("ape" %in% rownames(installed.packages()))) {
  install.packages("ape")
}
if (!("tidyverse" %in% rownames(installed.packages()))) {
  install.packages("tidyverse")
}

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
