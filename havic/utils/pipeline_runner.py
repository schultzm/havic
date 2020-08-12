#!/usr/bin/env python3

"""
Run the pipeline.

Go
"""

from .. import __ref_seq__, __parent_dir__
import pkg_resources
import io
import re
import sys
import os
from pathlib import Path, PurePath
from Bio import SeqIO
from ruffus import (mkdir,
                    follows,
                    files,
                    transform,
                    formatter,
                    jobs_limit,
                    pipeline_run,
                    pipeline_printout_graph,
                    originate)


def make_path(outdir, prefix, filename):
    """
    Make a filepath
    :param outdir:
    :param prefix:
    :param filename:
    :return: filepath
    """
    return Path(outdir).joinpath(f"{prefix}{filename}").as_posix()

class Pipeline:
    def __init__(self,
                 query_files,
                 trim_seqs,
                 subject_file,
                 redo,
                 n_snps,
                 seqlen,
                 matrixplots,
                 prefix,
                 outdir,
                 minimap2_kmer,
                 path_to_clusterpicker,
                 iqtree_threads):
        """
        Receive the arguments from argparse:

        :param query_files:
        :param trim_seqs:
        :param subject_file:
        :param redo:
        :param n_snps:
        :param seqlen:
        :param matrixplots:
        :param prefix:
        :param outdir:
        :param minimap2_kmer:
        """
        self.query_files = [Path(filename).resolve(strict=True) for filename in
                            query_files]
        self.trim_seqs = [re.sub('[^A-Za-z0-9]+', '_', i.replace("_(reversed)", "") \
                          .replace("(", "").replace(")", "").rstrip()) for i in trim_seqs]
        self.subject = subject_file
        if subject_file:
            self.subject = Path(self.subject).resolve(strict=True)
        else:
            self.subject = pkg_resources.resource_filename(__parent_dir__,
                                                           __ref_seq__)
        print(f"Will map amplicons to {self.subject}")
        self.refseq = SeqIO.read(self.subject, "fasta")
        self.reflen = len(self.refseq.seq)
        self.header = self.refseq.id
        self.redo = redo
        self.n_snps = n_snps
        self.seqlen = seqlen
        if matrixplots:
            self.matrixplots = "as.logical(TRUE)"
        else:
            self.matrixplots = "as.logical(FALSE)"
        self.prefix = prefix
        self.outdir = outdir
        self.outfiles = {
            'tmp_fasta': make_path(self.outdir, self.prefix,
                                   "HAV_all_tmpfasta.fa"),
            'tmp_bam': make_path(self.outdir,
                                 self.prefix, "HAV_all_minimap2.bam"),
            'tmp_bam_idx': make_path(self.outdir, self.prefix,
                                     "HAV_all_minimap2.bam.bai"),
            'bam2fasta': make_path(self.outdir, self.prefix,
                                   f"HAV_all_minimap2.bam2fasta.R"),
            'bam2fasta_Rout': make_path(self.outdir, self.prefix,
                                   f"HAV_all_minimap2.bam2fasta.Rout"),
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
                                            f".figTree"),
            'clusterpicked_tree_bkp': make_path(self.outdir, self.prefix,
                                                f"HAV_all_minimap2.stack"
                                                f".trimmed.fa.div_"
                                                f"{self.n_snps}SNPsIn"
                                                f"{self.seqlen}"
                                                f"bp.mp_clusterPicks"
                                                f".nwk.figTree"),
            'cluster_assignments': make_path(self.outdir, self.prefix,
                                             f"HAV_all_minimap2.stack.trimmed"
                                             f".fa.mp_clusterPicks_log.txt"),
            'clusters_assigned': make_path(self.outdir, self.prefix,
                                           f"HAV_all_minimap2.stack.trimmed"
                                           f".fa.mp_clusterPicks_summarised"
                                           f".txt"),
            'treeplotr': make_path(self.outdir, self.prefix,
                                   f"HAV_all_minimap2.stack.trimmed.fa"
                                   f".Rplot.R"),
            'treeplotr_out': make_path(self.outdir, self.prefix,
                                   f"HAV_all_minimap2.stack.trimmed.fa"
                                   f".Rplot.Rout")

        }
        self.minimap2_kmer = minimap2_kmer
        self.iqtree_threads = iqtree_threads
        self.path_to_clusterpicker = path_to_clusterpicker
        from .. import __ref_amplicon__
        self.havnet_ampliconseq = SeqIO.read(open(pkg_resources. \
                                           resource_filename(__parent_dir__,
                                                             __ref_amplicon__),
                                                  "r"), "fasta")

    def _compile_input_fasta(self):
        # 1 Compile the fasta files to single file
        quality_controlled_seqs = []
        # 1.01 Append the reference amplicon
        quality_controlled_seqs.append(self.havnet_ampliconseq)
        for query_file in self.query_files:
            print(query_file)
            for record in SeqIO.parse(query_file, "fasta"):
                record.id = re.sub('[^A-Za-z0-9]+', '_', record.id.replace("_(reversed)", "") \
                    .replace("(", "").replace(")", "").rstrip())
                # 1.02 Remove duplicates.
                if record.id not in [
                    record.id for record in quality_controlled_seqs
                ]:
                    quality_controlled_seqs.append(record)
                else:
                    print(f"Duplicate record found (only one copy of this "
                          f"added to quality_controlled_seqs): {record.id}")
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
        try:
            with open(self.outfiles['bam2fasta'], 'w') as out_r:
                from ..mapping.bam2fasta import bam2fasta
                cmd = bam2fasta % (
                    self.outfiles['tmp_bam'], self.outfiles['tmp_bam_idx'],
                    self.header, 1,
                    self.reflen,
                    self.outfiles['fasta_from_bam'])
                out_r.write(cmd)
            print(cmd)
            os.system(f"R CMD BATCH {self.outfiles['bam2fasta']} {self.outfiles['bam2fasta_Rout']}")
        except OSError:
            sys.exit("bam2fasta error.  Run 'havic check'.")

    def _get_clean_fasta_alignment(self):
        from Bio import AlignIO
        from Bio.Alphabet import generic_dna
        alignment = AlignIO.read(open(self.outfiles['fasta_from_bam'], "r"),
                                 "fasta", alphabet=generic_dna)
        # 4.1 Trim the alignment for isolates in arg.trim_seq to match
        # refamplicon.
        from ..utils.trim_alignment import Trimmed_alignment
        if not self.trim_seqs:
            self.trim_seqs = ''
        aln_trim = Trimmed_alignment(alignment,
                                     self.havnet_ampliconseq.id,
                                     '-',
                                     self.trim_seqs)
        if len(aln_trim.alignment) > 2:
            aln_trim._get_refseq_boundary()
            aln_trim.trim_seqs_to_ref()
            # 4.1.1 Depad the alignment.
            aln_trim.depad_alignment()
            AlignIO.write(aln_trim.alignment,
                          self.outfiles['fasta_from_bam_trimmed'],
                          "fasta")
        else:
            return aln_trim.alignment


    def _run_iqtree(self):
        if self.redo:
            redo = "redo -redo"
        else:
            redo = " TEST"
        cmd = f"iqtree -s {self.outfiles['fasta_from_bam_trimmed']} " \
              f"-T {self.iqtree_threads} -bb 1000 -m{redo}"
        # print(cmd)
        os.system(cmd)

    def _midpoint_root_iqtree(self):
        # 5.1 Midpoint root the phylogeny using ete3
        from ete3 import Tree
        tree = Tree(self.outfiles['treefile'], format=0)
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
        cmd = f"{self.path_to_clusterpicker} " \
              f"{self.outfiles['fasta_from_bam_trimmed']} " \
              f"{self.outfiles['mp_treefile']} 70 95 " \
              f"{self.n_snps/self.seqlen} 15 valid"
        print(cmd)
        os.system(cmd)

    def _bkp_clusterpickedtree(self):
        cmd = f"cp {self.outfiles['clusterpicked_tree']} " \
              f"{self.outfiles['clusterpicked_tree_bkp']}"
        print(cmd)
        os.system(cmd)

    def _summarise_cluster_assignments(self):
        from subprocess import Popen, PIPE
        import shlex
        cmd = f"grep ClusterNumber {self.outfiles['cluster_assignments']} -n"
        proc = Popen(shlex.split(cmd), stdout=PIPE)
        line_number = int(proc.communicate()[0].decode('UTF-8').split(":")[0])
        import pandas as pd
        df = pd.read_table(self.outfiles["cluster_assignments"],
                           skiprows=line_number - 1, header=0)
        with open(self.outfiles["clusters_assigned"], "w") as output_handle:
            output_handle.write(f"Isolate\tClusterNumber\n")
            for i in df.index.values:
                if isinstance(df.loc[i, "TipNames"], str):
                    tips = [j.replace("[", "").replace("]", "").strip() for j
                            in df.loc[i, "TipNames"].split(",")]
                    for tip in tips:
                        output_handle.write(
                            f"{tip}\tCluster_{df.loc[i, 'ClusterNumber']}\n")

    def _plot_results(self):
        """
        Link the alignment to the tree and plot it.

        :return: None
        """
        print("Starting results summaries using R")
        with open(self.outfiles['treeplotr'], "w") as out_r:
            from ..plotters.treeplot_snpplot import plot_functions
            cmd = plot_functions.replace("<- z", "<- \"" +
                                         self.outfiles[
                                             'fasta_from_bam_trimmed'] +
                                         "\"") \
                .replace("<- a", "<- " + str(self.n_snps)) \
                .replace("<- b", "<- " + str(self.seqlen)) \
                .replace("<- k", "<- " + str(self.minimap2_kmer)) \
                .replace("<- e", "<- " + str(self.matrixplots))
            # print(cmd)
            out_r.write(cmd)
        os.system(f"R CMD BATCH {self.outfiles['treeplotr']} {self.outfiles['treeplotr_out']}")

    def _run(self):
        """
        Run the pipeline using Ruffus.

        :return: None
        """

        # Pipeline starts here without Ruffus
        # For development and testing.
        # os.mkdir(self.outdir)
        # self._compile_input_fasta()
        # self._minimap2_input_fasta_to_ref()
        # self._bam2fasta()
        # self._get_clean_fasta_alignment()
        # self._run_iqtree()
        # self._midpoint_root_iqtree()
        # self._clusterpick()
        # self._plot_results()

        # if not os.path.exists(self.path_to_clusterpicker) or 'cluster' not in \
        #         self.path_to_clusterpicker.lower():
        #     sys.exit(f"ClusterPicker error: "
        #              f"Check {self.path_to_clusterpicker} exists and re-try.")

        # Pipeline starts here with Ruffus
        @mkdir(self.outdir)
        def create_outdir():
            pass
            # print(f"Creating output directory {self.outdir}")

        @follows(create_outdir)
        @files(self.query_files,
               self.outfiles['tmp_fasta'], self.havnet_ampliconseq)
        def compile_input_fasta(infile, outfile, refamplicon):
            self._compile_input_fasta()

        @follows(compile_input_fasta)
        @files(self.outfiles['tmp_fasta'],
               self.outfiles['tmp_bam'])
        def minimap2_input_fasta_to_ref(infile, outfile):
            self._minimap2_input_fasta_to_ref()

        @follows(minimap2_input_fasta_to_ref)
        @files(self.outfiles['tmp_bam'],
               self.outfiles['fasta_from_bam'])
        def bam2fasta(infile, outfile):
            self._bam2fasta()

        @follows(bam2fasta)
        @files(self.outfiles['fasta_from_bam'],
               self.outfiles['fasta_from_bam_trimmed'])
        def get_cleaned_fasta(infile, outfile):
            aln = self._get_clean_fasta_alignment()
            if aln and len(aln) < 3:
                exit_statement = f'{aln}\n' +\
                                 f'Need at least three sequences in ' +\
                                 f'alignment to continue (n={len(aln)})'
                import sys
                sys.exit(exit_statement)
        @follows(get_cleaned_fasta)
        @files(self.outfiles['fasta_from_bam_trimmed'],
               self.outfiles["mp_treefile"])
        def run_iqtree(infile, outfile):
            self._run_iqtree()

        @follows(run_iqtree)
        @files(self.outfiles['treefile'],
               self.outfiles['mp_treefile'])
        def midpoint_root_iqtree(infile, outfile):
            self._midpoint_root_iqtree()

        @follows(midpoint_root_iqtree)
        @files([self.outfiles['fasta_from_bam_trimmed'],
                self.outfiles['mp_treefile']],
               self.outfiles["clusterpicked_tree"])
        def clusterpick_from_mpr_iqtree_and_cleaned_fasta(infile, outfile):
            self._clusterpick()

        @follows(clusterpick_from_mpr_iqtree_and_cleaned_fasta)
        @files(self.outfiles["cluster_assignments"],
               self.outfiles["clusters_assigned"])
        def summarise_cluster_assignments(infile, outfile):
            self._summarise_cluster_assignments()

        @follows(clusterpick_from_mpr_iqtree_and_cleaned_fasta)
        @files(self.outfiles['clusterpicked_tree'],
               self.outfiles['clusterpicked_tree_bkp'])
        def backup_clusterpicked_figtree(infile, outfile):
            self._bkp_clusterpickedtree()

        @follows(clusterpick_from_mpr_iqtree_and_cleaned_fasta)
        @files([self.outfiles['fasta_from_bam_trimmed'],
                self.outfiles['mp_treefile']],
               self.outfiles['treeplotr'])
        def plot_results_ggtree(infiles, outfiles):
            self._plot_results()

        # Run the pipeline
        import tempfile

        def mv_tmp_sqlite(temp_sqlite_db, perm_sqlite_db):
            try:
                os.popen(f'mv {temp_sqlite_db} {perm_sqlite_db}')
                print(f"Moved {temp_sqlite_db} to {perm_sqlite_db}.")
            except IOError:
                print(f'Unable to move {temp_sqlite_db} to {self.outdir}',
                      file=sys.stderr)

        with tempfile.TemporaryDirectory() as tmpfile:
            db_name = '.ruffus_history.sqlite'
            temp_sqlite_db = Path(tmpfile).joinpath(db_name)
            perm_sqlite_db = Path(self.outdir).joinpath(db_name)
            if perm_sqlite_db.exists():
                os.popen(f'cp {perm_sqlite_db} {temp_sqlite_db}')
                print(f'Copied {perm_sqlite_db} to {temp_sqlite_db}.')
            else:
                print(f'Making new SQLite db at {temp_sqlite_db}')
            if self.redo:
                pipeline_run(forcedtorun_tasks=compile_input_fasta,
                             history_file=temp_sqlite_db)
                mv_tmp_sqlite(temp_sqlite_db, perm_sqlite_db)
            else:
                pipeline_run(history_file=temp_sqlite_db)
                mv_tmp_sqlite(temp_sqlite_db, perm_sqlite_db)

            # Print out the pipeline graph
            pipeline_printout_graph(
                make_path(self.outdir, self.prefix, "_pipeline_graph.svg"),
                "svg")
        # todo - 1.1 trim the sequences to remove primers
    #
    # def pipeline_of_pipelines(self):
    #     """
    #     Run the pipeline repeatedly.
    #     :return: None
    #     """
    #
    #     @originate(self.outdir)
    #     def run(self.outdir):
    #         self._run()
    #
    #     def run():
    #         self._run()
    #
    #     def run():
    #         self._run()
    #
    #     def run():
    #         self._run()
    #
    #     def run():
    #         self._run()
    #
    #     def run():
    #         self._run()
    #
    #     def run():
    #         self._run()
    #
    #     def run():
    #         self._run()
    #
    #     def run():
    #         self._run()
    #
    #     def run():
    #         self._run()
    #
    #     pipeline_run()


if __name__ == "__main__":
    import doctest

    doctest.testmod()
