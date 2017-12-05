bam2fasta = '''
library(GenomicAlignments)
library(Rsamtools)
bamFile <- BamFile('%s', index = '%s')
bamFile
seqinfo(bamFile)
outp <- stackStringsFromBam(bamFile, param=GRanges('%s', IRanges(%d, %d)), Lpadding.letter='-', Rpadding.letter='-', use.names = T)
sink('%s')
for(i in 1:length(outp)){
    cat('>',names(outp[i]), '\n', sep='')
    cat(toString(outp[[i]]), '\n', sep='')
}
sink()
'''