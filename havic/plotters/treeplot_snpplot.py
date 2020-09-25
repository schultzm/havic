# todo Need to allow for labelling 'new' tips in tree figure
plot_functions = """
library(phytools)
library(tidyverse)
library(ggtree)
library(Biostrings)
library(colorspace)

basename <- z
print(paste0('basename: ', basename))
distfract <- a
supportvals <- b
method <- hh
kmer <- k
matrixplots <- e
highlight <- wz
tree <- read.newick(file=paste0(basename, '.rooted.treefile'))
tree_clust <- read.newick(file=paste0(basename, '.rooted_clusterPicks.nwk'))
cluster_picks <- matrix(nrow = 0, ncol = 2)
colnames(cluster_picks) <- c('Isolate', 'Cluster')
cluster_picks

for(tip in tree_clust$tip.label[grep('Clust', tree_clust$tip.label)]){
    outp_list_of_lists <- strsplit(tip, split = '_')
    output_list <- outp_list_of_lists[[1]]
    cluster_picks <- rbind(cluster_picks,
                           c(paste0(output_list[2:length(output_list)],
                           collapse='_'), output_list[1]))
    }

cluster_picks <- data.frame(cluster_picks, stringsAsFactors=FALSE)
cluster_picks


list_of_clusters <- split(cluster_picks$Isolate, cluster_picks$Cluster)

if(matrixplots){
    plt <- ggtree(tree, size=0.1) %<+% cluster_picks
    plt$data <- plt$data %>% add_column(Highlight = as.character(NA))
    for(i in highlight) {
        plt$data[which(plt$data$label == i), 'Highlight'] <- 'query'
    }
    plt$data$Highlight <- factor(plt$data$Highlight, levels=c('context', 'query'))
    plt$data$labclust <- as.character(NA)
    for(row in 1:nrow(plt$data)) {
        if(is.na(plt$data[row, 'Cluster'])) {
            plt$data[row, 'labclust'] <- plt$data[row, 'label']
        }
        else {
            plt$data[row, 'labclust'] <- paste0(plt$data[row, 'label'], '_', plt$data[row, 'Cluster'])
        }
    }
    print(as.data.frame(plt$data))
    str(plt$data)
    offst <- max(dist.nodes(tree))/4
    # if(offst <= 0){
    #     offst <- 0.1
    # }
    fntsz <- -0.005093*length(tree$tip.label)+1.955556
    if(fntsz < 0.5) {
        fntsz <- 0.5
        }
    print(paste('FontSize for plot', fntsz))

    wdth <- 1
    q <- plt + geom_tiplab(aes(label=label,
                               color=Cluster),
                           align=TRUE,
                           size=fntsz,
                           linesize=0.1) +
        geom_tippoint(aes(shape=Highlight), size=fntsz, color='red', alpha=0.7) +
        geom_text2(aes(x=branch,
                       label=as.integer(label),
                       vjust=-0.3,
                       hjust=1,
                       subset=(isTip!=TRUE & as.integer(label)>=70),
                       na.rm=TRUE),
                   size=fntsz,
                   na.rm=TRUE) +
        geom_treescale(fontsize = fntsz,
                       linesize=0.1) +
        ggtitle(label = "ML IQtree with bootstrap %, tips cluster-picked (left); fasta alignment (right)",
                subtitle = paste0('Clusters (coloured labels) have been picked as clades with >= ',
                                  supportvals,
                                  '% support and divergence <= ',
                                  distfract*100,
                                  '% (distance method=\\'',method, '\\')')
               ) +
        scale_color_discrete_qualitative(palette = "Dark3", na.value = "black") +
        guides(color = guide_legend(override.aes = list(size = 3)))
    # q
    pdf(file=paste0(basename,
                   '.rooted.treefile_',
                   distfract*100,
                   'percent_divergence_',
                   method, '_msa.pdf'),
        paper = 'a4r',
        width=11.69,
        height=8.27)

    h <- msaplot(p=q,
                 fasta=basename,
                 offset = offst,
                 width = wdth,
                 bg_line = FALSE,
                 color=c('#f7fcfd', #white
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
                         '#f8c0ff')
                )
    print(h)
    dev.off()
}


aln <- read.dna(basename, format = 'fasta', as.character = TRUE)
snp_dists <- function(alignment, exclude_char){
    mat <- matrix(0, nrow = nrow(alignment), ncol = nrow(alignment))
    colnames(mat) <- rownames(alignment)
    rownames(mat) <- rownames(alignment)
    mat_snps <- mat
    pw_dist <- function(Aln_sub){
        len_aln <- Aln_sub[, c(which(!(Aln_sub[1,] %in% exclude_char) &
                                     !(Aln_sub[2,] %in% exclude_char)
                                    )
                              ),
                           drop=FALSE]
        n_snps <- len_aln[, c(which(len_aln[1,]!=len_aln[2,])), drop=FALSE]
        return(c(paste0('=', ncol(n_snps), '/', ncol(len_aln)), ncol(n_snps)))
    }
    for(i in 1:nrow(alignment)){
        for(j in i:nrow(alignment)){
            seq1 <- rownames(alignment)[i]
            seq2 <- rownames(alignment)[j]
            val <- pw_dist(alignment[c(i,j), , drop=FALSE])
            mat[seq1, seq2] <- val[1]
            mat[seq2, seq1] <- val[1]
            mat_snps[seq1, seq2] = as.numeric(val[2])
            mat_snps[seq2, seq1] = as.numeric(val[2])
            }
        }
    return(list(mat, mat_snps))
}

exclusions <- tolower(c(names(IUPAC_CODE_MAP)[!(names(IUPAC_CODE_MAP) %in%
                                                c('A', 'C', 'T', 'G')
                                               )
                                             ],
                        '-', '?')
                     )

#Returns two matrices
snps <- snp_dists(aln, exclusions)
heatmap_data <- snps[[2]]

if(matrixplots){
    library(pheatmap)
    pdf(file = paste0(basename,
                      '_SNPdists.pdf'),
        width=11.69,
        height=8.27)
    annos <- cluster_picks
    rownames(annos) <- cluster_picks$Isolate
    annos <- annos['Cluster']
    labs_row <- c()
    for(i in 1:length(row.names(heatmap_data))){
        labs_row <- c(labs_row, paste0(annos[row.names(heatmap_data)[i],'Cluster'], '_', row.names(heatmap_data)[i]))
    }
    colors_anfang <- unique(ggplot_build(q)$data[[3]][,c('colour', 'label')])
    # colors_anfang$cluster <- as.character(NA)
    colnames(colors_anfang) <- c('Colour', 'Isolate')
    cluster_colors <- unique(merge(cluster_picks, colors_anfang, by='Isolate')[,c('Cluster', 'Colour')])
    cluster_colors_ <- cluster_colors$Colour
    names(cluster_colors_) <- cluster_colors$Cluster
    cluster_colors_ <- list(Cluster=cluster_colors_)
    print(pheatmap(heatmap_data,
         show_rownames = TRUE,
         show_colnames = TRUE,
         annotation_names_row = TRUE,
         annotation_names_col = TRUE,
         annotation_row = annos,
         annotation_col = annos,
         annotation_colors = cluster_colors_,
         labels_row = labs_row,
         display_numbers = TRUE,
         number_format = '%.0f',
         fontsize = fntsz*2))
    dev.off()
}


write.csv(x=heatmap_data, file=paste0(basename,
                                      '_SNPdists.csv'),
          quote=FALSE)

write.csv(snps[[1]], file=paste0(basename,
                                 '_SNPcountsOverAlignLength.csv'),
          quote = FALSE)
"""
