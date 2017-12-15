ggtree_plot = '''
library(phytools)
library(plotly)
library(tidyverse)
library(ggtree)

basename <- z
print(paste0('basename', basename))
nsnps <- a
seqlen <- b
kmer <- k
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

fntsz <- -0.00458*length(tree$tip.label)+2.23817
wdth <- 4
q <- p + geom_tiplab(aes(label=paste0(label, ' ', Cluster), color=Cluster), size=fntsz, linesize=0.2, align=TRUE) +
    theme(legend.position = "right") +
    geom_tippoint(aes(color=Cluster), size=fntsz, na.rm=T) +
    geom_text2(aes(x=branch, label=as.integer(label), vjust=-0.3, hjust=1, subset=(isTip!=TRUE), na.rm=TRUE), size=fntsz, na.rm=TRUE) +
    geom_treescale(x=0.01, y =-2, offset=1, fontsize = fntsz) +
    annotate("text", x = 0.015, y=-4, label = "Substitutions per site", size=fntsz) + 
    ggtitle(label = "ML IQtree with bootstrap %, tips cluster-picked (left); fasta alignment (right)", subtitle = paste0('Clusters have been picked at 95% support allowing ', nsnps, ' SNPs in ', seqlen, ' bp'))

pdf(file=paste0(basename, '.mp.treefile.', 'div_', nsnps, 'SNPsIn', seqlen, 'bp_msa.pdf'), paper = 'a4r', width=11.69, height=8.27, onefile = TRUE)

h <- msaplot(p=q, fasta=basename, offset = 0.2, width = wdth, bg_line = FALSE, color=c('#f7fcfd', #white
                                                                                    '#ef3b2c', #red
                                                                                    '#41ab5d', #green
                                                                                    '#ffffbf', #yellow
                                                                                    '#4292c6', #blue
                                                                                    '#dface5',
                                                                                    '#c699cc',
                                                                                    '#ad86b2',
                                                                                    '#947399',
                                                                                    '#7c607f',
                                                                                    '#634c66',
                                                                                    '#4a394c',
                                                                                    '#312633',
                                                                                    '#181319',
                                                                                    '#000000',
                                                                                    '#f8c0ff'))
h
dev.off()


library(pheatmap)
heatmap_data <- read.FASTA(basename)
heatmap_data <- dist.dna(heatmap_data, model='N', as.matrix = TRUE)
re_name <- function(iname, trmatrix) {
    if(iname %in% trmatrix[,1]){
        row <- which(trmatrix[,1]==iname)
        return(paste0(trmatrix[row,1], ' ', trmatrix[row,2])) 
    }
    else{
        return(iname)
    }
}

colnames(heatmap_data) <- lapply(colnames(heatmap_data),
                                 re_name,
                                 cluster_picks)
rownames(heatmap_data) <- lapply(rownames(heatmap_data),
                                 re_name,
                                 cluster_picks)

pheatmap(heatmap_data,
         show_rownames = T,
         show_colnames = T,
         display_numbers = TRUE,
         number_format = '%.0f',
         fontsize = fntsz*2,
         filename = paste0(basename,
                           '_ClusterNamesPicked',
                           nsnps,
                           'SNPsIn',
                           seqlen,
                           '_SNPdists.pdf'),
         width=11.69,
         height=8.27)
dev.off()
write.csv(x=heatmap_data, file=paste0(basename,
                                      '_ClusterNamesPicked',
                                      nsnps,
                                      'SNPsIn',
                                      seqlen,
                                      '_SNPdists.csv'),
          quote=FALSE)

'''