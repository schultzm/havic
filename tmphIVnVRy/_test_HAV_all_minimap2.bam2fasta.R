
library(GenomicAlignments)
library(Rsamtools)
bamFile <- BamFile('/Users/mschultz/Desktop/MDU_DAMG/MDU/HAVNET/HAVIC/tmphIVnVRy/_test_HAV_all_minimap2.bam', index = '/Users/mschultz/Desktop/MDU_DAMG/MDU/HAVNET/HAVIC/tmphIVnVRy/_test_HAV_all_minimap2.bam.bai')
bamFile
seqinfo(bamFile)
outp <- stackStringsFromBam(bamFile,
                            param=GRanges('NC_001489.1', IRanges(1, 7478)),
                            Lpadding.letter='-',
                            Rpadding.letter='-',
                            use.names = T)
sink('/Users/mschultz/Desktop/MDU_DAMG/MDU/HAVNET/HAVIC/tmphIVnVRy/_test_HAV_all_minimap2.stack.fa')
for(i in 1:length(outp)){
    cat('>',names(outp[i]), '
', sep='')
    cat(toString(outp[[i]]), '
', sep='')
}
sink()
