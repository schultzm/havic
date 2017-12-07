plottree = '''
library(phytools)

tree <- read.newick('/Users/mschultz/HAV_all_minimap2.stack.trimmed.fa.treefile')
outgrp <- NULL
rt_tree <- function(tree, outgroup){
    if(is.null(outgrp)){
        tree <- ladderize(midpoint.root(tree), right = FALSE)
    }
    else{
        tree <-root(tree, outgroup=outgroup,resolve.root=T)
    }
    return(tree)
}

tree <- rt_tree(tree, outgrp)

library(tidyverse)
library(ggtree)


ggtree_tre <- ggtree(tree) + geom_tiplab(size=0.9, align=TRUE, linesize = 0.2) + geom_text2(aes(subset=!isTip, label=label, vjust=-1, hjust=2), size=0.9)
pdf('/Users/mschultz/HAV_all_minimap2.stack.trimmed.fa.treefile.msa.pdf', paper = 'a4r', width=11.69, height=8.27)
msaplot(p=ggtree_tre, fasta='/Users/mschultz/HAV_all_minimap2.stack.trimmed.fa', offset = 0.1, width = 11.69, bg_line = FALSE)
dev.off()
''''