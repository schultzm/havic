ggtree_plot = '''
library(phytools)
library(plotly)
library(tidyverse)
library(ggtree)

basename <- z
print(paste0('basename', basename))
divergence <- w
tree <- read.newick(file=paste0(basename, '.mp.treefile'))
tree_clust <- read.newick(file=paste0(basename, '.mp_clusterPicks.nwk'))
cluster_picks <- matrix(nrow = 0, ncol = 2) 
colnames(cluster_picks) <- c('Isolate', 'Cluster')
cluster_picks

for(tip in tree_clust$tip.label[grep('Clust', tree_clust$tip.label)]){
    outp_list_of_lists <- strsplit(tip, split = '_') 
    output_list <- outp_list_of_lists[[1]]
    cluster_picks <- rbind(cluster_picks, c(paste0(output_list[2:length(output_list)], collapse='_'), output_list[1]))
    }

cluster_picks <- data.frame(cluster_picks)
cluster_picks

list_of_clusters <- split(cluster_picks$Isolate, cluster_picks$Cluster)
names(list_of_clusters)
p <- ggtree(tree) %<+% cluster_picks
names(p)
q <- p + geom_tiplab(aes(label=label, color=Cluster), size=1, linesize=0.2, align=TRUE) +
    theme(legend.position = "right") +
    geom_tippoint(aes(color=Cluster), size=0.5, na.rm=T) +
    geom_text2(aes(x=branch, label=as.integer(label), vjust=-0.3, hjust=1, subset=(isTip!=TRUE), na.rm=TRUE), size=1, na.rm=TRUE) +
    geom_treescale(x=0.01, y =-2, offset=1, fontsize = 1) +
    annotate("text", x = 0.015, y=-4, label = "Substitutions per site", size=1) + ggtitle(label = "ML IQtree with bootstrap %, tips cluster-picked (left); fasta alignment (right)", subtitle = paste0('Clusters picked at ', divergence, ' divergence'))
h <- msaplot(p=q, fasta=basename, offset = 0.1, width = 11.69, bg_line = FALSE)
pdf(file=paste0(basename, '.mp.treefile.', 'div_', divergence, '_msa.pdf'), paper = 'a4r', width=11.69, height=8.27, onefile = TRUE)
h
dev.off()

'''