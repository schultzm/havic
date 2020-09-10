#!/usr/bin/env python3

"""
Run the pipeline.

Go
"""

from .. import __parent_dir__  # , __nextflow_nf__
from pkg_resources import resource_filename as rf
import re
import sys
import os
from pathlib import Path
import pandas as pd
import shutil
from subprocess import Popen, PIPE
import shlex
from Bio import SeqIO
from ruffus import (
    mkdir,
    follows,
    files,
    pipeline_run,
    pipeline_printout_graph as pipeprintgraph,
)


def make_path(parentdir, filename):
    """
    Make a filepath
    :param parentdir: a parent directory
    :param filename: a file in the parent dir
    :return: joined filepath
    """
    return Path(parentdir).joinpath(filename).as_posix()


def absolute_path(fname_in, default_path):
    """Get absolute paths for filename.

    Args:
        fname (string): filename
        test_status(boolean): This file is a pre-packaged havic datafile.

    Returns:
        valid absolute file path if it is a file, else None
    """
    fname_out = None
    if default_path:
        fname_in = rf(__parent_dir__, fname_in)
    else:
        pass
    if Path(fname_in).is_file():
        fname_out = Path(fname_in).resolve(strict=True)
    else:
        print(f"\nWarning, '{fname_in}' is not a valid file path.\n")
    return fname_out


def correct_characters(input_string):
    """Remove non alphanumeric characters from string.

    Args:
        input_string (string)

    Returns:
        string: output string with non alpha-numeric characters removed.
    """
    output_string = re.sub(
        "[^A-Za-z0-9]+",
        "_",
        input_string.replace("_(reversed)", "")
        .replace("(", "")
        .replace(")", "")
        .rstrip(),
    )
    return output_string


class Pipeline:
    def __init__(self, yaml_in):
        """Read the dictionary, and make it available to Pipeline() methods.
        The dictionary contains all the run parameters, such as minimap2
        or iqtree2 settings.

        Args:
            yaml_in (dict): A dictionary object parsed from a yaml input file.
        """
        self.yaml_in = yaml_in  # a dict object
        self.query_files = list(
            filter(
                None,
                [
                    absolute_path(fname, yaml_in["DEFAULT_QUERIES"])
                    for fname in yaml_in["QUERY_FILES"]
                ],
            )
        )
        if not self.query_files:
            sys.exit("Unable to continue without input query_files.")
        self.trim_seqs = [correct_characters(i) for i in yaml_in["TRIM_SEQS"]]
        self.subject = absolute_path(yaml_in["SUBJECT_FILE"], yaml_in["DEFAULT_REFS"])
        for key, value in yaml_in.items():
            print(key, value)
        self.refseq = SeqIO.read(self.subject, "fasta")
        self.reflen = len(self.refseq.seq)
        self.header = self.refseq.id
        self.outdir = yaml_in["OUTDIR"]
        repstr = yaml_in["RUN_PREFIX"]
        self.outfiles = {
            "tmp_fasta": make_path(self.outdir, f"{repstr}tmpfasta.fa"),
            "seq_header_replacements": make_path(
                self.outdir, f"{repstr}seq_id_replace.tsv"
            ),
            "duplicates": make_path(self.outdir, f"{repstr}duplicate_seqs.txt"),
            "tmp_bam": make_path(self.outdir, f"{repstr}minimap2.bam"),
            "tmp_bam_idx": make_path(self.outdir, f"{repstr}minimap2.bam.bai"),
            "bam2fasta": make_path(self.outdir, f"{repstr}minimap2.bam2fasta.R"),
            "bam2fasta_Rout": make_path(
                self.outdir, f"{repstr}minimap2.bam2fasta.Rout"
            ),
            "fasta_from_bam": make_path(self.outdir, f"{repstr}minimap2.stack.fa"),
            "fasta_from_bam_trimmed": make_path(
                self.outdir, f"{repstr}minimap2.stack.trimmed.fa"
            ),
            "treefile": make_path(
                self.outdir, f"{repstr}minimap2.stack.trimmed.fa.treefile"
            ),
            "mp_treefile": make_path(
                self.outdir, f"{repstr}minimap2.stack.trimmed.fa.mp.treefile"
            ),
            "clusterpicked_tree": make_path(
                self.outdir,
                f"{repstr}minimap2.stack.trimmed.fa.mp_clusterPicks.nwk.figTree",
            ),
            "clusterpicked_tree_bkp": make_path(
                self.outdir,
                f"{repstr}minimap2.stack"
                f".trimmed.fa.div_"
                f"{yaml_in['CLUSTER_PICKER_SETTINGS']['distance_fraction']}"
                f"distancefraction"
                f".mp_clusterPicks"
                f".nwk.figTree",
            ),
            "cluster_assignments": make_path(
                self.outdir,
                f"{repstr}minimap2.stack.trimmed.fa.mp_clusterPicks_log.txt",
            ),
            "clusters_assigned": make_path(
                self.outdir,
                f"{repstr}minimap2.stack.trimmed"
                f".fa.mp_clusterPicks_summarised"
                f".txt",
            ),
            "treeplotr": make_path(
                self.outdir, f"{repstr}minimap2.stack.trimmed.fa.Rplot.R"
            ),
            "treeplotr_out": make_path(
                self.outdir, f"{repstr}minimap2.stack.trimmed.fa.Rplot.Rout"
            ),
        }

        self.minimap2_cmd = (
            f"minimap2 -k {yaml_in['MINIMAP2_SETTINGS']['k_mer']} "
            f"-a {self.subject} "
            f"{self.outfiles['tmp_fasta']} "
            f"| samtools sort > {self.outfiles['tmp_bam']}"
        )
        self.clusterpick_cmd = (
            f"{yaml_in['CLUSTER_PICKER_SETTINGS']['executable']} "
            f"{self.outfiles['fasta_from_bam_trimmed']} "
            f"{self.outfiles['mp_treefile']} "
            f"{yaml_in['CLUSTER_PICKER_SETTINGS']['coarse_subtree_support']} "
            f"{yaml_in['CLUSTER_PICKER_SETTINGS']['fine_cluster_support']} "
            f"{yaml_in['CLUSTER_PICKER_SETTINGS']['distance_fraction']} "
            f"{yaml_in['CLUSTER_PICKER_SETTINGS']['large_cluster_threshold']} "
            f"{yaml_in['CLUSTER_PICKER_SETTINGS']['distance_method']}"
        )
        self.iqtree_cmd = str(
            f"{yaml_in['IQTREE2_SETTINGS']['executable']} "
            f"-s {self.outfiles['fasta_from_bam_trimmed']} "
            f"{yaml_in['IQTREE2_SETTINGS']['threads']} "
            f"{yaml_in['IQTREE2_SETTINGS']['model_finder']}"
            f"{yaml_in['IQTREE2_SETTINGS']['state_frequency']} "
            f"{yaml_in['IQTREE2_SETTINGS']['ultrafast_bootstrap']} "
            f"{yaml_in['IQTREE2_SETTINGS']['protect_violations']} "
            f"{yaml_in['IQTREE2_SETTINGS']['redo']}"
        )
        self.target_region = SeqIO.read(
            open(absolute_path(yaml_in["SUBJECT_TARGET_REGION"], yaml_in["DEFAULT_REFS"]), "r"), "fasta"
        )

    def _compile_input_fasta(self):
        # 1 Compile the fasta files to single file
        quality_controlled_seqs = []
        # 1.01 Append the reference amplicon
        quality_controlled_seqs.append(self.target_region)
        keyval_ids = {}
        dups = []
        for query_file in self.query_files:
            for record in SeqIO.parse(query_file, "fasta"):
                input_id = record.id
                record.id = correct_characters(record.id)
                # if str(record.id) != str(input_id):
                keyval_ids[str(input_id)] = str(record.id)
                # 1.02 Remove duplicates.
                if record.id not in [record.id for record in quality_controlled_seqs]:
                    quality_controlled_seqs.append(record)
                else:
                    dups.append(str(record.id))
        if keyval_ids:
            with open(self.outfiles["seq_header_replacements"], "w") as out_h:
                out_h.write("INPUT_SEQ_HEADER\tOUTPUT_SEQ_HEADER\n")
                for key, val in keyval_ids.items():
                    out_h.write(key + "\t" + val + "\n")
        else:
            print("Zero seq headers were modified.")
        if dups:
            with open(self.outfiles["duplicates"], "w") as out_h:
                out_h.write("\n".join(dups))
        else:
            print("Zero duplicate sequences were found.")
        SeqIO.write(quality_controlled_seqs, self.outfiles["tmp_fasta"], "fasta")

    def _minimap2_input_fasta_to_ref(self):
        cmd = self.minimap2_cmd
        print(cmd)
        os.system(cmd)
        cmd = f"samtools index {self.outfiles['tmp_bam']}"
        os.system(cmd)
        # 3.1 find the unmapped sequences.
        cmd = f"samtools view -f 4 {self.outfiles['tmp_bam']}"
        cmd2 = "cut -f 1"
        proc = Popen(shlex.split(cmd), stdout=PIPE, stderr=PIPE)
        proc2 = Popen(shlex.split(cmd2), stdin=proc.stdout, stdout=PIPE, stderr=PIPE)
        result = proc2.communicate()[0].decode("UTF-8").split("\n")
        if result:
            print(
                f"\nUnmapped reads at k-mer "
                f"{self.yaml_in['MINIMAP2_SETTINGS']['k_mer']}:"
            )
            print("\n".join(result))
        else:
            pass

    def _bam2fasta(self):
        """
        Convert the bam file to fasta by stacking strings on ref to get MSA.
        :return: MSA fasta from input bam file
        """
        with open(self.outfiles["bam2fasta"], "w") as out_r:
            from ..mapping.bam2fasta import bam2fasta

            cmd = bam2fasta % (
                self.outfiles["tmp_bam"],
                self.outfiles["tmp_bam_idx"],
                self.header,
                1,
                self.reflen,
                self.outfiles["fasta_from_bam"],
            )
            out_r.write(cmd)
        # print(cmd)
        os.system(f"R CMD BATCH {self.outfiles['bam2fasta']} {self.outfiles['bam2fasta_Rout']}")

    def _get_clean_fasta_alignment(self):
        from Bio import AlignIO
        from Bio.Alphabet import generic_dna

        alignment = AlignIO.read(
            open(self.outfiles["fasta_from_bam"], "r"), "fasta", alphabet=generic_dna
        )
        # 4.1 Trim the alignment for isolates in arg.trim_seq to match
        # refamplicon.
        from ..utils.trim_alignment import Trimmed_alignment

        if not self.trim_seqs:
            self.trim_seqs = ""
        aln_trim = Trimmed_alignment(
            alignment, self.target_region.id, "-", self.trim_seqs
        )
        if len(aln_trim.alignment) > 2:
            aln_trim._get_refseq_boundary()
            aln_trim.trim_seqs_to_ref()
            # 4.1.1 Depad the alignment.
            aln_trim.depad_alignment()
            AlignIO.write(
                aln_trim.alignment, self.outfiles["fasta_from_bam_trimmed"], "fasta"
            )
        else:
            return aln_trim.alignment

    def _run_iqtree(self):
        os.system(self.iqtree_cmd)

    def _midpoint_root_iqtree(self):
        # 5.1 Midpoint root the phylogeny using ete3
        from ete3 import Tree

        tree = Tree(self.outfiles["treefile"], format=0)
        root = tree.get_midpoint_outgroup()
        tree.set_outgroup(root)
        tree.ladderize(direction=1)
        # dist_formatter is to prevent scientific notation.
        # with branch lengths in scientific notation, ClusterPicker dies.
        tree.write(outfile=self.outfiles["mp_treefile"], dist_formatter="%0.16f")

    def _clusterpick(self):
        """
        Run CLUSTER_PICKER on the tree and alignment
        :return: None
        """
        os.system(self.clusterpick_cmd)

    def _bkp_clusterpickedtree(self):
        shutil.copyfile(
            self.outfiles["clusterpicked_tree"], self.outfiles["clusterpicked_tree_bkp"]
        )

    def _summarise_cluster_assignments(self):
        cmd = f"grep ClusterNumber {self.outfiles['cluster_assignments']} -n"
        proc = Popen(shlex.split(cmd), stdout=PIPE)
        line_number = int(proc.communicate()[0].decode("UTF-8").split(":")[0])
        df = pd.read_table(
            self.outfiles["cluster_assignments"], skiprows=line_number - 1, header=0
        )
        with open(self.outfiles["clusters_assigned"], "w") as output_handle:
            output_handle.write("Isolate\tClusterNumber\n")
            for i in df.index.values:
                if isinstance(df.loc[i, "TipNames"], str):
                    tips = [
                        j.replace("[", "").replace("]", "").strip()
                        for j in df.loc[i, "TipNames"].split(",")
                    ]
                    for tip in tips:
                        output_handle.write(
                            f"{tip}\tCluster_{df.loc[i, 'ClusterNumber']}\n"
                        )

    def _plot_results(self):
        """
        Link the alignment to the tree and plot it.

        :return: None
        """
        print("Starting results summaries using R")
        with open(self.outfiles["treeplotr"], "w") as out_r:
            from ..plotters.treeplot_snpplot import plot_functions
            tiphighlights = "c('" + "', '".join([correct_characters(i) for i in self.yaml_in["HIGHLIGHT_TIP"]]) + "')"
            cmd = (
                plot_functions.replace(
                    "basename <- z",
                    'basename <- "' + self.outfiles["fasta_from_bam_trimmed"] + '"',
                )
                .replace(
                    "distfract <- a",
                    "distfract <- "
                    + str(self.yaml_in["CLUSTER_PICKER_SETTINGS"]["distance_fraction"]),
                )
                .replace(
                    "supportvals <- b",
                    "supportvals <- "
                    + str(
                        self.yaml_in["CLUSTER_PICKER_SETTINGS"]["fine_cluster_support"]
                    ),
                )
                .replace(
                    "method <- hh",
                    'method <- "'
                    + str(self.yaml_in["CLUSTER_PICKER_SETTINGS"]["distance_method"])
                    + '"',
                )
                .replace(
                    "kmer <- k",
                    "kmer <- " + str(self.yaml_in["MINIMAP2_SETTINGS"]["k_mer"]),
                )
                .replace(
                    "matrixplots <- e",
                    "matrixplots <- " + str(self.yaml_in["PLOTS"]).upper(),
                )
                .replace(
                    "highlight <- wz",
                    "highlight <- " + tiphighlights,
                )
            )
            # print(cmd)
            out_r.write(cmd)
        os.system(f"R CMD BATCH {self.outfiles['treeplotr']} {self.outfiles['treeplotr_out']}")

    def _run(self):
        """
        Run the pipeline using Ruffus.

        :return: None
        """

        # os.system(f"nextflow run {rf(__parent_dir__, __nextflow_nf__)}")
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

        # Pipeline starts here with Ruffus
        @mkdir(self.outdir)
        def create_outdir():
            pass

        @follows(create_outdir)
        @files(self.query_files, self.outfiles["tmp_fasta"], self.target_region)
        def compile_input_fasta(infile, outfile, refamplicon):
            self._compile_input_fasta()

        @follows(compile_input_fasta)
        @files(self.outfiles["tmp_fasta"], self.outfiles["tmp_bam"])
        def minimap2_input_fasta_to_ref(infile, outfile):
            self._minimap2_input_fasta_to_ref()

        @follows(minimap2_input_fasta_to_ref)
        @files(self.outfiles["tmp_bam"], self.outfiles["fasta_from_bam"])
        def bam2fasta(infile, outfile):
            self._bam2fasta()

        @follows(bam2fasta)
        @files(self.outfiles["fasta_from_bam"], self.outfiles["fasta_from_bam_trimmed"])
        def get_cleaned_fasta(infile, outfile):
            aln = self._get_clean_fasta_alignment()
            if aln and len(aln) < 3:
                exit_statement = (
                    f"{aln}\n"
                    f"Need at least three sequences in "
                    f"alignment to continue (n={len(aln)})"
                )
                import sys

                sys.exit(exit_statement)

        @follows(get_cleaned_fasta)
        @files(self.outfiles["fasta_from_bam_trimmed"], self.outfiles["mp_treefile"])
        def run_iqtree(infile, outfile):
            self._run_iqtree()

        @follows(run_iqtree)
        @files(self.outfiles["treefile"], self.outfiles["mp_treefile"])
        def midpoint_root_iqtree(infile, outfile):
            self._midpoint_root_iqtree()

        @follows(midpoint_root_iqtree)
        @files(
            [self.outfiles["fasta_from_bam_trimmed"], self.outfiles["mp_treefile"]],
            self.outfiles["clusterpicked_tree"],
        )
        def clusterpick_from_mpr_iqtree_and_cleaned_fasta(infile, outfile):
            self._clusterpick()

        @follows(clusterpick_from_mpr_iqtree_and_cleaned_fasta)
        @files(self.outfiles["cluster_assignments"], self.outfiles["clusters_assigned"])
        def summarise_cluster_assignments(infile, outfile):
            self._summarise_cluster_assignments()

        @follows(clusterpick_from_mpr_iqtree_and_cleaned_fasta)
        @files(
            self.outfiles["clusterpicked_tree"], self.outfiles["clusterpicked_tree_bkp"]
        )
        def backup_clusterpicked_figtree(infile, outfile):
            self._bkp_clusterpickedtree()

        @follows(clusterpick_from_mpr_iqtree_and_cleaned_fasta)
        @files(
            [self.outfiles["fasta_from_bam_trimmed"], self.outfiles["mp_treefile"]],
            self.outfiles["treeplotr"],
        )
        def plot_results_ggtree(infiles, outfiles):
            self._plot_results()

        # Run the pipeline
        import tempfile

        with tempfile.TemporaryDirectory() as tmpfile:
            db_name = ".ruffus_history.sqlite"
            temp_sqlite = Path(tmpfile).joinpath(db_name)
            perm_sqlite = Path(self.outdir).joinpath(db_name)
            if self.yaml_in["FORCE_OVERWRITE_AND_RE_RUN"]:
                if Path(self.outdir).exists():
                    shutil.rmtree(Path(self.outdir))
                else:
                    pass
                pipeline_run(forcedtorun_tasks=create_outdir, history_file=temp_sqlite)
                shutil.copyfile(temp_sqlite, perm_sqlite)
            else:
                shutil.copyfile(perm_sqlite, temp_sqlite)
                pipeline_run(history_file=temp_sqlite)
                shutil.copyfile(temp_sqlite, perm_sqlite)

            # Print out the pipeline graph
            pipeprintgraph(make_path(self.outdir, "pipeline_graph.svg"), "svg")

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
