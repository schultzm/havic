#!/usr/bin/env python3

# minimap2 -x asm20 -k 5 NC_001489.fa example.fa --MD -a | samtools sort | samtools view -h > test.bam
# samtools index test.bam

import pysam
import sys
from collections import defaultdict

file = sys.argv[1]
alnfile = pysam.AlignmentFile(file, 'rb')
from Bio import SeqIO, AlignIO
from Bio.AlignIO import MultipleSeqAlignment as MSA
from Bio.Alphabet import generic_dna
# print((alnfile.lengths))
# print(help(alnfile.pileup))
aln = defaultdict(list)

print((alnfile.lengths[0]))
sys.exit()
for aligned_segment in alnfile.fetch():
    cigarred_seq = []
    aln_dict = aligned_segment.to_dict()
    # minimap2 stores a reverse complemented seq in sam output
    aln_pairs = aligned_segment.get_aligned_pairs(with_seq=True)
    start_pad, start_iter, end_pad = False, 0, 0
    while start_pad is False:
        for query_pos, ref_pos, ref_base in aln_pairs:
            if isinstance(ref_pos, type(None)):
                pass
            else:
                start_pad, start_iter = ref_pos, query_pos
    cigarred_seq.append(start_pad*'-')
    for query_pos, ref_pos, ref_nt in aln_pairs[start_iter:]:
        if isinstance(ref_pos, type(None)):
            pass
        elif isinstance(query_pos, type(None)):
            cigarred_seq.append('-')
        else:
            cigarred_seq.append(aln_dict['seq'][query_pos])
            if query_pos < alnfile.lengths[0]
    print(f">{aln_dict['name']}\n{len(''.join(cigarred_seq))}")
