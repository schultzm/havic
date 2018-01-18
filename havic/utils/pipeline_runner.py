#!/usr/bin/env python3

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
from Bio import SeqIO
from ruffus import (mkdir,
                    follows,
                    files,
                    transform,
                    formatter,
                    jobs_limit,
                    pipeline_run,
                    pipeline_printout_graph)


def make_path(outdir, prefix, filename):
    """
    Make a filepath
    :param outdir:
    :param prefix:
    :param filename:
    :return: filepath
    """
    import os
    return os.path.join(os.path.abspath(outdir), prefix + filename)


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
        self.query_files = [Input_file(file, "Query").filename for file in
                            query_files]
        self.trim_seqs = trim_seqs

        if subject_file is None:
            self.subject = pkg_resources.resource_filename(__parent_dir__,
                                                           __ref_seq__)
        else:
            self.subject = Input_file(self.subject_file, "Subject").filename
        print(f"Will map amplicons to {self.subject}")

        self.refseq = SeqIO.read(self.subject, "fasta")
        self.reflen = len(self.refseq.seq)
        self.header = self.refseq.id
        self.redo = redo
        self.n_snps = n_snps
        self.seqlen = seqlen
        self.prefix = prefix
        self.outdir = outdir
        self.outfiles = {
            'tmp_fasta': make_path(self.outdir, self.prefix,
                                   "HAV_all_tmpfasta.fa"),
            'tmp_bam': make_path(self.outdir,
                                 self.prefix, "HAV_all_minimap2.bam"),
            'tmp_bam_idx': make_path(self.outdir, self.prefix,
                                     "HAV_all_minimap2.bam.bai"),
            'fasta_from_bam': make_path(self.outdir, self.prefix,
                                        "HAV_all_minimap2.stack.fa"),
            'fasta_from_bam_trimmed': make_path(self.outdir, self.prefix,
                                                "HAV_all_minimap2.stack"
                                                f".trimmed.fa"),
            'treefile': make_path(self.outdir, self.prefix,
                                  f"HAV_all_minimap2.stack.trimmed.fa"
                                  f".treefile"),
            'mp_treefile': make_path(self.outdir, self.prefix,
                                     f"HAV_all_minimap2.stack.trimmed"
                                     f".fa.mp.treefile"),
            'clusterpicked_tree': make_path(self.outdir, self.prefix,
                                            f"HAV_all_minimap2.stack.trimmed"
                                            f".fa.mp_clusterPicks.nwk"
                                            f".FigTree"),
            'clusterpicked_tree_bkp': make_path(self.outdir, self.prefix,
                                                f"HAV_all_minimap2.stack"
                                                f".trimmed.fa.div_"
                                                f"{self.n_snps}SNPsIn"
                                                f"{self.seqlen}"
                                                f"bp.mp_clusterPicks"
                                                f".nwk.figTree"),
            'treeplotr': make_path(self.outdir, self.prefix,
                                   f"HAV_all_minimap2.stack.trimmed.fa"
                                   f".Rplot.R")

        }

        self.minimap2_kmer = minimap2_kmer
        from ..data.havnet_amplicon import havnet_ampliconseq
        self.havnet_ampliconseq = havnet_ampliconseq

    def _compile_input_fasta(self):
        # 1 Compile the fasta files to single file
        import os
        from Bio import SeqIO
        quality_controlled_seqs = []
        # if not os.path.exists(os.path.abspath(self.outdir)):
        #     os.mkdir(os.path.abspath(self.outdir))
        for query_file in self.query_files:
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
        quality_controlled_seqs.append(
            SeqIO.read(io.StringIO(self.havnet_ampliconseq), "fasta"))
        SeqIO.write(quality_controlled_seqs, self.outfiles['tmp_fasta'],
                    "fasta")

    def _minimap2_input_fasta_to_ref(self):
        cmd = f"minimap2 -k {self.minimap2_kmer} -a {self.subject} " \
              f"{self.outfiles['tmp_fasta']} " \
              f"| samtools sort > {self.outfiles['tmp_bam']}"
        print(cmd)
        os.system(cmd)
        cmd = f"samtools index {self.outfiles['tmp_bam']}"
        os.system(cmd)
        # 3.1 find the unmapped sequences.
        cmd = f"samtools view -f 4 {self.outfiles['tmp_bam']} | cut -f 1"
        print(f"Unmapped reads at k-mer {self.minimap2_kmer}:")
        os.system(cmd)

    def _bam2fasta(self):
        """

        :return: fasta from input bam file
        """
        from rpy2 import robjects
        import warnings
        from rpy2.rinterface import RRuntimeWarning
        warnings.filterwarnings("ignore", category=RRuntimeWarning)
        try:
            from ..mapping.bam2fasta import bam2fasta
            cmd = bam2fasta % (
                self.outfiles['tmp_bam'], f"{self.outfiles['tmp_bam']}.bai",
                self.header, 1,
                self.reflen,
                self.outfiles['fasta_from_bam'])
            print(cmd)
            robjects.r(cmd)
        except OSError:
            sys.exit("bam2fasta error")

    def _get_clean_fasta_alignment(self):
        from Bio import AlignIO
        alignment = AlignIO.read(open(self.outfiles['fasta_from_bam'], "r"),
                                 "fasta")
        # 4.1 Trim the alignment for isolates in arg.trim_seq to match
        # refamplicon.
        from ..utils.trim_alignment import Trimmed_alignment
        aln_trim = Trimmed_alignment(alignment,
                                     SeqIO.read(
                                         io.StringIO(self.havnet_ampliconseq),
                                         "fasta").id, '-', self.trim_seqs)
        aln_trim._get_refseq_boundary()
        aln_trim.trim_seqs_to_ref()
        # 4.1.1 Depad the alignment.
        aln_trim.depad_alignment()
        AlignIO.write(aln_trim.alignment,
                      self.outfiles['fasta_from_bam_trimmed'],
                      "fasta")

    def _run_iqtree(self):
        if self.redo:
            redo = "redo -redo"
        else:
            redo = " TEST"
        cmd = f"iqtree -s {self.outfiles['fasta_from_bam_trimmed']} " \
              f"-nt AUTO -bb 1000 -m{redo}"
        print(cmd)
        os.system(cmd)

    def _midpoint_root_iqtree(self):
        # 5.1 Midpoint root the phylogeny using ete3
        from ete3 import Tree
        print(f"Reading {self.outfiles['treefile']}")
        tree = Tree(self.outfiles['treefile'], format=0)
        print(tree.write())
        root = tree.get_midpoint_outgroup()
        tree.set_outgroup(root)
        tree.ladderize(direction=1)
        # tree.resolve_polytomy(default_dist=0.01)
        print(tree.write())
        # dist_formatter is to prevent scientific notation.
        # with branch lengths in scientific notation, ClusterPicker dies.
        tree.write(outfile=self.outfiles["mp_treefile"],
                   dist_formatter="%0.16f")

    def _clusterpick(self):
        """
        Run CLUSTER_PICKER on the tree and alignment
        :return: None
        """
        from ..tests.dependencies import CLUSTER_PICKER
        cmd = f"java -jar {CLUSTER_PICKER} " \
              f"{self.outfiles['fasta_from_bam_trimmed']} " \
              f"{self.outfiles['mp_treefile']} 70 95 " \
              f"{self.n_snps/self.seqlen} 15 valid"
        print(cmd)
        os.system(cmd)
        cmd = f"cp {self.outfiles['clusterpicked_tree']} " \
              f"{self.outfiles['clusterpicked_tree_bkp']}"
        os.system(cmd)

    def _plot_results(self):
        """
        Link the alignment to the tree and plot it.

        :return: None
        """
        with open(f"{self.outfiles['fasta_from_bam_trimmed']}.Rplot.R",
                  "w") as out_r:
            from ..plotters.treeplot_snpplot import plot_functions
            cmd = plot_functions.replace("<- z", "<- \"" +
                                         self.outfiles[
                                             'fasta_from_bam_trimmed'] +
                                         "\"") \
                .replace("<- a", "<- " + str(self.n_snps)) \
                .replace("<- b", "<- " + str(self.seqlen)) \
                .replace("<- k", "<- " + str(self.minimap2_kmer))
            print(cmd)
            out_r.write(cmd)
        os.system(f"Rscript {self.outfiles['fasta_from_bam_trimmed']}.Rplot.R")

    def run(self):
        """
        Run the pipeline using Ruffus.

        :return: None
        """

        # Pipeline starts here
        @mkdir(self.outdir)
        def create_outdir():
            print(f"Creating output directory {self.outdir}")

        @follows(create_outdir)
        @transform(self.query_files, formatter(None),
                   self.outfiles['tmp_fasta'], self.havnet_ampliconseq)
        def compile_input_fasta(infile, outfile, refamplicon):
            self._compile_input_fasta()

        #
        @follows(compile_input_fasta)
        @transform(self.outfiles['tmp_fasta'], formatter(None),
                   self.outfiles['tmp_bam'])
        def minimap2_input_fasta_to_ref(infile, outfile):
            self._minimap2_input_fasta_to_ref()

        @follows(minimap2_input_fasta_to_ref)
        @transform(self.outfiles['tmp_bam'], formatter(None),
                   self.outfiles['fasta_from_bam'])
        def bam2fasta(infile, outfile):
            self._bam2fasta()

        @follows(bam2fasta)
        @transform(self.outfiles['fasta_from_bam'], formatter(None),
                   self.outfiles['fasta_from_bam_trimmed'])
        def get_cleaned_fasta(infile, outfile):
            self._get_clean_fasta_alignment()

        @follows(get_cleaned_fasta)
        @transform(self.outfiles['fasta_from_bam_trimmed'], formatter(None),
                   self.outfiles["mp_treefile"])
        def run_iqtree(infile, outfile):
            self._run_iqtree()

        @follows(run_iqtree)
        @transform(self.outfiles['treefile'], formatter(None),
                   self.outfiles['mp_treefile'])
        def midpoint_root_iqtree(infile, outfile):
            self._midpoint_root_iqtree()

        @follows(midpoint_root_iqtree)
        @transform([self.outfiles['fasta_from_bam_trimmed'],
                    self.outfiles['mp_treefile']],
                   formatter(None),
                   self.outfiles["clusterpicked_tree"])
        def clusterpick_from_iqtree_and_cleaned_fasta(infiles, outfile):
            self._clusterpick()

        @follows(clusterpick_from_iqtree_and_cleaned_fasta)
        @transform([self.outfiles['fasta_from_bam_trimmed'],
                    self.outfiles['mp_treefile'],
                    self.outfiles['clusterpicked_tree']],
                   formatter(None), self.outfiles['treeplotr'])
        def plot_results(infiles, outfiles):
            self._plot_results()

        # Run the pipeline
        if self.redo:
            pipeline_run(forcedtorun_tasks=compile_input_fasta)
        else:
            pipeline_run()

        # Print out the pipeline graph
        pipeline_printout_graph(
            make_path(self.outdir, self.prefix, "_pipeline_graph.svg"),
            "svg")
        # todo - 1.1 trim the sequences to remove primers


if __name__ == "__main__":
    import doctest
    doctest.testmod()
