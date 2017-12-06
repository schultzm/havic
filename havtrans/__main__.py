#!/usr/bin/env python3

'''
Parse fasta files (remove spaces in descriptor etc)
Concatenate to single file
Run the pipeline
Trim the alignment
'''

from havtrans.classes.classes import Input_file, Check_dependency
from havtrans.tests import check_r_dependencies
from havtrans.tests.dependencies import SOFTWAREZ, R_LIBS, CLUSTER_PICKER
from havtrans.mapping.bam2fasta import bam2fasta
import sys

def main():
    '''
    The main routine.
    '''
    import argparse
    import os
    import sys
    parser = argparse.ArgumentParser(description='Run HAVTrans')
    subparsers = parser.add_subparsers(title='Sub-commands help',
                                       help='', metavar='',
                                       dest='subparser_name')
    subparser_run = subparsers.add_parser('run', help='Run the analysis.',
                                           description='Run the pipeline.')
    subparser_run.add_argument('-q', '--query_files', help='Query file',
                               nargs='+', required=True)
    subparser_run.add_argument('-s', '--subject_file', help='Subject file',
                               required=True)
    subparser_version = subparsers.add_parser('version', help='Print version.',
                                           description='Print version.')
    args = parser.parse_args()
    
    if not args.subparser_name:
        os.system('havtrans -h')

    queries = [Input_file(file, 'Query').filename for file in args.query_files]
    subject = Input_file(args.subject_file, 'Subject').filename

#     for dep in SOFTWAREZ:
#         path = Check_dependency(dep)
#         path.check_software()
#     for dep in R_LIBS: #move this to class
#         check_r_dependencies.importr_tryhard(dep)
#         print(f'R library {dep}'.ljust(28)+': ok', file=sys.stderr)

    
    #1 Compile the query fasta files to single file
    from Bio import SeqIO
    quality_controlled_seqs = []
    tmp_fasta = os.path.expanduser('~/all_tmp.fa')
    tmp_bam = os.path.expanduser('~/HAV_all_minimap2.bam')
    fasta_from_bam = os.path.expanduser('~/HAV_all_minimap2.stack.fa')
    fasta_from_bam_trimmed = os.path.expanduser('~/HAV_all_minimap2.stack.trimmed.fa')
    for query_file in queries:
        for record in SeqIO.parse(query_file, 'fasta'):
            record.id = record.id.replace('_(reversed)', ' reversed') \
                              .replace('(', '').replace(')', '')
            quality_controlled_seqs.append(record)
    SeqIO.write(quality_controlled_seqs, tmp_fasta, 'fasta')
    
    #1.1 trim the sequences to remove primers

    #2 get ref and ref stats
    refseq = SeqIO.read(subject, 'fasta')
    reflen = len(refseq.seq)
    header = refseq.id
    
    #3 get minimap2 done
    import os
    cmd = f'minimap2 -k 28 -a {subject} {tmp_fasta} | samtools sort > {tmp_bam}'
    cmd2 = f'samtools index {tmp_bam}'
    os.system(cmd)
    os.system(cmd2)


    #4 get the fasta from the bam using bam2fasta
    from rpy2 import robjects
    import warnings
    from rpy2.rinterface import RRuntimeWarning
    warnings.filterwarnings("ignore", category=RRuntimeWarning)
    try:
        robjects.r(bam2fasta % (tmp_bam, f'{tmp_bam}.bai', header, 1, reflen, fasta_from_bam))
    except:
        sys.exit('bam2fasta error')
    
    #4.1 Trim the alignment
    from Bio import AlignIO
    alignment = AlignIO.read(open(fasta_from_bam, 'r'), 'fasta')
    site_set = {'-'}
    start_pos = 0
    while len(site_set) == 1:
        site_set = set(alignment[:, start_pos])
        start_pos += 1
    site_set = {'-'}
    end_pos = alignment.get_alignment_length()-1 #subtract one due to 0 and 1-based indexing issues
    while len(site_set) == 1:
        site_set = set(alignment[:, end_pos])
        end_pos -= 1
    alignment_trimmed = alignment[:, start_pos:end_pos+2]#.format('fasta') #Add 2, again due to indexing discrepancies
    AlignIO.write(alignment_trimmed, fasta_from_bam_trimmed, 'fasta')
    #5 Run iqtree on the extracted bam2fasta
    cmd = f'iqtree -s {fasta_from_bam_trimmed} -nt AUTO -bb 1000 -m TEST'
    os.system(cmd)
    
    #6 Run CLUSTER_PICKER on the tree and alignment
    cmd = f'java -jar {CLUSTER_PICKER} {fasta_from_bam_trimmed} {fasta_from_bam_trimmed}.treefile 70 95 0.006 15 valid'
    print(cmd)
    os.system(cmd)
    
    #6 Link tree to alignment
    from ete3 import PhyloTree, TreeStyle
    #Load a tree and link it to an alignment.
    t = PhyloTree(open(f'{fasta_from_bam_trimmed}.treefile', 'r').read())
    t.ladderize()
    ts = TreeStyle()
    ts.show_branch_support = True
    ts.show_branch_length = True
    t.link_to_alignment(alignment=alignment_trimmed.format('fasta'), alg_format='fasta')
    t.render(f'{fasta_from_bam_trimmed}.treefile.pdf',
             tree_style=ts)

if __name__ == '__main__':
    main()

