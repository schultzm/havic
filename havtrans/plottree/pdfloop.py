looper = '''
library(phytools)
library(ape)
library(tidyverse)
library(ggtree)

tree1 <- read.newick('%s')
tree <- read.newick('%s')

cluster_picks <- matrix(nrow = 0, ncol = 2) 
colnames(cluster_picks) <- c('Cluster', 'Isolate')
cluster_picks

for(tip in tree$tip.label[grep('Clust', tree$tip.label)]){
    outp_list_of_lists <- strsplit(tip, split = '_') 
    output_list <- outp_list_of_lists[[1]]
    cluster_picks <- rbind(cluster_picks, c(output_list[1], paste0(output_list[2:length(output_list)], collapse='_')))
    }

cluster_picks <- data.frame(cluster_picks, stringsAsFactors=FALSE)
cluster_picks

list_of_clusters <- split(cluster_picks, cluster_picks$Cluster)
names(list_of_clusters)


plot_list <- list()
for(cluster_n in 1:length(list_of_clusters)){
    cluster <- as.matrix(list_of_clusters[cluster_n][[1]])
    print(cluster)
    nrow(cluster)
        tips = c(cluster[,'Isolate'])
        print(tips)
        nde <- findMRCA(tree1, tips, type = 'node')
        
        tre <- ape::extract.clade(tree1, node=nde)
        print(nde)
        titl <- paste(names(list_of_clusters[cluster_n]), ': branch_lengths=\'none\' (cladogram)')
        ggtree_tre <- ggtree(tre, branch.length = 'none') +
              geom_tiplab(size=2, align=TRUE, linesize = 0.2) +
              geom_text2(aes(subset=!isTip, label=label, vjust=-1, hjust=2), size=2) +
              ggtitle(titl)

    
        q <- msaplot(p=ggtree_tre, fasta='/Users/mschultz/HAV_all_minimap2.stack.trimmed.fa', bg_line = FALSE, width = 100, offset = 10) 
        plot_list[[cluster_n]] <- q
    }

dev.off()
pdf(file = paste('/Users/mschultz/HAV_all_minimap2.stack.trimmed.fa.treefile', '.msa.pdf', sep=''), paper = 'a4r', width=11.69, height=8.27, onefile = TRUE)
for(i in plot_list){
    print(i)
}
dev.off()
'''