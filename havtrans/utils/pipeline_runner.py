"""
Run the pipeline.

Go
"""

from .. import __ref_seq__, __parent_dir__
from ..utils.input_file import Input_file
import pkg_resources
import io
import sys
import os
from ruffus import (mkdir,
                    follows,
                    files,
                    jobs_limit,
                    pipeline_run,
                    pipeline_printout_graph, )


class Pipeline:
    def __init__(self,
                 query_files,
                 trim_seqs,
                 subject_file,
                 redo,
                 n_snps,
                 seqlen,
                 prefix,
                 outdir,
                 minimap2_kmer):
        """
        Receive the arguments from argparse:

        :param query_files:
        :param trim_seqs:
        :param subject_file:
        :param redo:
        :param n_snps:
        :param seqlen:
        :param prefix:
        :param outdir:
        :param minimap2_kmer:
        """
        self.query_files = query_files
        self.trim_seqs = trim_seqs
        self.subject_file = subject_file
        self.redo = redo
        self.n_snps = n_snps
        self.seqlen = seqlen
        self.prefix = prefix
        self.outdir = outdir
        self.minimap2_kmer = minimap2_kmer

    def piperun(self):
        """

        :return: None
        """
        import os
        self.queries = [Input_file(file, "Query").filename for file in
                   self.query_files]
        print('trim seqs', self.trim_seqs, 'queries', self.queries)
        if self.subject_file is None:
            subject = pkg_resources.resource_filename(__parent_dir__,
                                                      __ref_seq__)
        else:
            subject = Input_file(self.subject_file, "Subject").filename
        print(f"Will map amplicons to {subject}")

        # 1 Compile the fasta files to single file
        from Bio import SeqIO
        quality_controlled_seqs = []
        if not os.path.exists(os.path.abspath(self.outdir)):
            os.mkdir(os.path.abspath(self.outdir))
        outfiles = {
            'tmp_fasta': os.path.join(os.path.abspath(self.outdir),
                                      self.prefix + "all_tmp.fa"),
            'tmp_bam': os.path.join(os.path.abspath(self.outdir),
                                    self.prefix + "HAV_all_minimap2.bam"),
            'fasta_from_bam': os.path.join(
                os.path.abspath(self.outdir),
                self.prefix + "HAV_all_minimap2.stack.fa"),
            'fasta_from_bam_trimmed': os.path.join(
                os.path.abspath(self.outdir),
                self.prefix + "HAV_all_minimap2" + ".stack.trimmed.fa")
        }
        for query_file in self.queries:
            print(query_file)
            for record in SeqIO.parse(query_file, "fasta"):
                record.id = record.id.replace("_(reversed)", "") \
                    .replace("(", "").replace(")", "")
                # 1.02 Remove duplicates.
                if record.id not in [
                    record.id for record in quality_controlled_seqs
                ]:
                    quality_controlled_seqs.append(record)
                else:
                    print(f"Duplicate record found (only one copy of this " +
                          "added to quality_controlled_seqs): {record.id}")
        # 1.01 Append the reference amplicon
        from ..data.havnet_amplicon import havnet_ampliconseq
        quality_controlled_seqs.append(
            SeqIO.read(io.StringIO(havnet_ampliconseq), "fasta"))
        SeqIO.write(quality_controlled_seqs, outfiles['tmp_fasta'], "fasta")
        # todo - 1.1 trim the sequences to remove primers
        # 2 get ref and ref stats
        refseq = SeqIO.read(subject, "fasta")
        reflen = len(refseq.seq)
        header = refseq.id
        # 3 get minimap2 done
        import os
        cmd = f"minimap2 -k {self.minimap2_kmer} -a {subject} " \
              f"{outfiles['tmp_fasta']} " \
              f"| samtools sort > {outfiles['tmp_bam']}"
        print(cmd)
        os.system(cmd)
        cmd = f"samtools index {outfiles['tmp_bam']}"
        os.system(cmd)
        # 3.1 find the unmapped sequences.
        cmd = f"samtools view -f 4 {outfiles['tmp_bam']} | cut -f 1"
        print(f"Unmapped reads at k-mer {self.minimap2_kmer}:")
        os.system(cmd)
        # 4 get the fasta from the bam using bam2fasta
        from rpy2 import robjects
        import warnings
        from rpy2.rinterface import RRuntimeWarning
        warnings.filterwarnings("ignore", category=RRuntimeWarning)
        try:
            from ..mapping.bam2fasta import bam2fasta
            cmd = bam2fasta % (
                outfiles['tmp_bam'], f"{outfiles['tmp_bam']}.bai", header, 1,
                reflen,
                outfiles['fasta_from_bam'])
            print(cmd)
            robjects.r(cmd)
        except OSError:
            sys.exit("bam2fasta error")
        # 4.0 Import the alignment
        from Bio import AlignIO
        alignment = AlignIO.read(open(outfiles['fasta_from_bam'], "r"),
                                 "fasta")
        # 4.1 Trim the alignment for isolates in arg.trim_seq to match
        # refamplicon.
        from ..utils.trim_alignment import Trimmed_alignment
        aln_trim = Trimmed_alignment(alignment,
                                     SeqIO.read(
                                         io.StringIO(havnet_ampliconseq),
                                         "fasta").id, '-', self.trim_seqs)
        aln_trim._get_refseq_boundary()
        aln_trim.trim_seqs_to_ref()
        # 4.1.1 Depad the alignment.
        aln_trim.depad_alignment()
        AlignIO.write(aln_trim.alignment, outfiles['fasta_from_bam_trimmed'],
                      "fasta")
        # 5 Run iqtree on the extracted bam2fasta
        if self.redo:
            redo = "redo -redo"
        else:
            redo = " TEST"
        cmd = f"iqtree -s {outfiles['fasta_from_bam_trimmed']} -nt AUTO -bb " \
              f"1000 -m{redo}"
        # TN+I+G4{redo}"
        os.system(cmd)
        # 5.1 Midpoint root the phylogeny using ete3
        treefile = f"{outfiles['fasta_from_bam_trimmed']}.treefile"

        print(f"Reading {treefile}")
        mp_treefile = f"{outfiles['fasta_from_bam_trimmed']}.mp.treefile"
        from ete3 import Tree
        tree = Tree(treefile, format=0)
        print(tree.write())
        root = tree.get_midpoint_outgroup()
        tree.set_outgroup(root)
        tree.ladderize(direction=1)
        # tree.resolve_polytomy(default_dist=0.01)
        print(tree.write())
        # dist_formatter is to prevent scientific notation.
        # with branch lengths in scientific notation, ClusterPicker dies.
        tree.write(outfile=mp_treefile, dist_formatter="%0.16f")
        # 6 Run CLUSTER_PICKER on the tree and alignment
        from ..tests.dependencies import CLUSTER_PICKER
        cmd = f"java -jar {CLUSTER_PICKER} " \
              f"{outfiles['fasta_from_bam_trimmed']} " \
              f"{mp_treefile} 70 95 {self.n_snps/self.seqlen} 15 valid"
        print(cmd)
        os.system(cmd)
        cmd = f"cp {outfiles['fasta_from_bam_trimmed']}.mp_clusterPicks.nwk" \
              f".figTree " \
              f"{outfiles['fasta_from_bam_trimmed']}.div_{self.n_snps}" \
              f"SNPsIn{self.seqlen}" \
              f"bp.mp_clusterPicks.nwk.figTree"
        os.system(cmd)
        # 6 Link tree to alignment and plot it
        # treestring = open(f"{outfiles[
        # 'fasta_from_bam_trimmed'].mp.treefile", "r").read()
        with open(f"{outfiles['fasta_from_bam_trimmed']}.Rplot.R",
                  "w") as out_r:
            from ..plotters.treeplot_snpplot import plot_functions
            cmd = plot_functions.replace("<- z", "<- \"" +
                                         outfiles[
                                             'fasta_from_bam_trimmed'] +
                                         "\"") \
                .replace("<- a", "<- " + str(self.n_snps)) \
                .replace("<- b", "<- " + str(self.seqlen)) \
                .replace("<- k", "<- " + str(self.minimap2_kmer))
            print(cmd)
            out_r.write(cmd)
        os.system(f"Rscript {outfiles['fasta_from_bam_trimmed']}.Rplot.R")
