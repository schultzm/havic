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
        .replace(":", "_")
        .rstrip(),
    )
    return output_string


class Pipeline:
    def __init__(self, yaml_in):
        """Read the dictionary, and make it available to Pipeline() methods.
        The dictionary contains all the run parameters, such as mapping
        or iqtree2 settings.

        Args:
            yaml_in (dict): A dictionary object parsed from a yaml input file.
        """
        self.yaml_in = yaml_in  # a dict object
        # yaml_in["TRIM_SEQS"] = [item for sublist in yaml_in["TRIM_SEQS"] for item in sublist]
        for key, value in yaml_in.items():
            print(f"{key}: {value}")
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
        self.trim_requests = {correct_characters(i): i for i in yaml_in["TRIM_SEQS"]} # this allows printing of unfound TRIM_SEQS
        self.trim_seqs = list(set(filter(None, [correct_characters(i) for i in yaml_in["TRIM_SEQS"]])))
        self.subject = absolute_path(yaml_in["SUBJECT_FILE"], yaml_in["DEFAULT_SUBJECT"])
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
            "tmp_bam": make_path(self.outdir, f"{repstr}map.bam"),
            "tmp_bam_idx": make_path(self.outdir, f"{repstr}map.bam.bai"),
            "bam2fasta": make_path(self.outdir, f"{repstr}map.bam2fasta.R"),
            "bam2fasta_Rout": make_path(
                self.outdir, f"{repstr}map.bam2fasta.Rout"
            ),
            "fasta_from_bam": make_path(self.outdir, f"{repstr}map.stack.fa"),
            "fasta_from_bam_trimmed": make_path(
                self.outdir, f"{repstr}map.stack.trimmed.fa"
            ),
            "treefile": make_path(
                self.outdir, f"{repstr}map.stack.trimmed.fa.treefile"
            ),
            "rooted_treefile": make_path(
                self.outdir, f"{repstr}map.stack.trimmed.fa.rooted.treefile"
            ),
            "clusterpicked_tree": make_path(
                self.outdir,
                f"{repstr}map.stack.trimmed.fa.rooted_clusterPicks.nwk.figTree",
            ),
            "cluster_assignments": make_path(
                self.outdir,
                f"{repstr}map.stack.trimmed.fa.rooted_clusterPicks_log.txt",
            ),
            "treeplotr": make_path(
                self.outdir, f"{repstr}map.stack.trimmed.fa.Rplot.R"
            ),
            "treeplotr_out": make_path(
                self.outdir, f"{repstr}map.stack.trimmed.fa.Rplot.Rout"
            ),
        }

        self.map_cmd = (
            f"{yaml_in['MAPPER_SETTINGS']['executable']} {yaml_in['MAPPER_SETTINGS']['other']} {yaml_in['MAPPER_SETTINGS']['k_mer']} "
            f"-a {self.subject} "
            f"{self.outfiles['tmp_fasta']} "
            f"| samtools view -h -F 256 -F 2048 | samtools sort > {self.outfiles['tmp_bam']}"
        )
        self.clusterpick_cmd = (
            f"{yaml_in['CLUSTER_PICKER_SETTINGS']['executable']} "
            f"{self.outfiles['fasta_from_bam_trimmed']} "
            f"{self.outfiles['rooted_treefile']} "
            f"{yaml_in['CLUSTER_PICKER_SETTINGS']['coarse_subtree_support']} "
            f"{yaml_in['CLUSTER_PICKER_SETTINGS']['fine_cluster_support']} "
            f"{yaml_in['CLUSTER_PICKER_SETTINGS']['distance_fraction']} "
            f"{yaml_in['CLUSTER_PICKER_SETTINGS']['large_cluster_threshold']} "
            f"{yaml_in['CLUSTER_PICKER_SETTINGS']['distance_method']}"
        )
        self.iqtree_cmd = str(
            f"{yaml_in['IQTREE2_SETTINGS']['executable']} "
            f"-s {self.outfiles['fasta_from_bam_trimmed']} "
            f"{yaml_in['IQTREE2_SETTINGS']['other']}"
        )
        self.target_region = SeqIO.read(
            open(absolute_path(yaml_in["SUBJECT_TARGET_REGION"], yaml_in["DEFAULT_SUBJECT"]), "r"), "fasta"
        )
        self.target_region.id = correct_characters(self.target_region.id)
        self.target_region.seq = self.target_region.seq.ungap("-")
        self.root = correct_characters(self.yaml_in["TREE_ROOT"])
        self.replacedheaders = None

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
                    record.seq = record.seq.ungap("-")
                    quality_controlled_seqs.append(record)
                else:
                    dups.append(str(record.id))
        if keyval_ids:
            self.replacedheaders = keyval_ids
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
        if self.root != 'midpoint' and self.root not in [record.id for record in quality_controlled_seqs]:
            lbreak = "\n"
            sys.exit(f"Incorrect specification of tree root (Hint: must be either 'midpoint' or sample from input fasta, but was {self.root}.  Choices are:{lbreak}{lbreak.join([record.id for record in quality_controlled_seqs])}")
        else:
            SeqIO.write(quality_controlled_seqs, self.outfiles["tmp_fasta"], "fasta")

    def _map_input_fasta_to_ref(self):
        cmd = self.map_cmd
        print(cmd)
        os.system(cmd)
        cmd = f"samtools index {self.outfiles['tmp_bam']}"
        os.system(cmd)
        # Find and print the unmapped sequences.
        cmd = f"samtools view -f 4 {self.outfiles['tmp_bam']}"
        cmd2 = "cut -f 1"
        proc = Popen(shlex.split(cmd), stdout=PIPE, stderr=PIPE)
        proc2 = Popen(shlex.split(cmd2), stdin=proc.stdout, stdout=PIPE, stderr=PIPE)
        result = proc2.communicate()[0].decode("UTF-8").split("\n")
        if result:
            print(
                f"\nUnmapped reads at k-mer "
                f"{self.yaml_in['MAPPER_SETTINGS']['k_mer']}:"
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
        """Give the alignment a haircut.

        Returns:
            MSA: The Biopython Multiple Sequence Alignment object
        """
        from Bio import AlignIO

        alignment = AlignIO.read(
            open(self.outfiles["fasta_from_bam"], "r"), "fasta")
        from ..utils.trim_alignment import Trimmed_alignment
        aln_trim = Trimmed_alignment(
            alignment, self.target_region.id, "-", self.trim_seqs
        )

        notfound = set(self.trim_seqs) - set([seq.id for seq in alignment])
        if notfound: # report trim requests in yaml_in['TRIM_SEQS'] that could not be found and trimmed
            for nf in notfound:
                print(f"Unable to find and trim {self.trim_requests[nf]}")

        if len(aln_trim.alignment) < 3:
            sys.exit('Not enough sequences to perform analysis.  Exiting now.\n')
        aln_trim.get_refseq_boundary()
        aln_trim.trim_seqs_to_ref()
        aln_trim.depad_alignment()
        AlignIO.write(
            aln_trim.alignment, self.outfiles["fasta_from_bam_trimmed"], "fasta"
        )
        return aln_trim.alignment

    def _run_iqtree(self):
        os.system(self.iqtree_cmd)

    def root_iqtree(self):
        """Midpoint or user-defined root setting of iqtree.
        """
        from ete3 import Tree
        tree = Tree(self.outfiles["treefile"], format=0)
        root_ = self.root
        root = None
        if root_ == 'midpoint':
            root = tree.get_midpoint_outgroup()
        else:
            root = root_
        tree.set_outgroup(root)
        tree.ladderize(direction=1)
        # dist_formatter is to prevent scientific notation.
        # with branch lengths in scientific notation, ClusterPicker dies.
        tree.write(outfile=self.outfiles["rooted_treefile"], dist_formatter="%0.16f")

    def _clusterpick(self):
        """
        Run CLUSTER_PICKER on the tree and alignment
        :return: None
        """
        os.system(self.clusterpick_cmd)

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
                    "kmer <- " + str(self.yaml_in["MAPPER_SETTINGS"]["k_mer"].rstrip().split(' ')[-1]),
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
        def map_input_fasta_to_ref(infile, outfile):
            self._map_input_fasta_to_ref()

        @follows(map_input_fasta_to_ref)
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
        @files(self.outfiles["fasta_from_bam_trimmed"], self.outfiles["rooted_treefile"])
        def run_iqtree(infile, outfile):
            self._run_iqtree()

        @follows(run_iqtree)
        @files(self.outfiles["treefile"], self.outfiles["rooted_treefile"])
        def root_iqtree(infile, outfile):
            self.root_iqtree()

        @follows(root_iqtree)
        @files(
            [self.outfiles["fasta_from_bam_trimmed"], self.outfiles["rooted_treefile"]],
            self.outfiles["clusterpicked_tree"],
        )
        def clusterpick_from_rooted_iqtree_and_cleaned_fasta(infile, outfile):
            self._clusterpick()

        @follows(clusterpick_from_rooted_iqtree_and_cleaned_fasta)
        @files(
            [self.outfiles["fasta_from_bam_trimmed"], self.outfiles["rooted_treefile"]],
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
                for fname in Path(self.outdir).glob(f"{self.yaml_in['RUN_PREFIX']}*"):
                    Path.unlink(fname)
                pipeline_run(forcedtorun_tasks=create_outdir, history_file=temp_sqlite)
                shutil.copyfile(temp_sqlite, perm_sqlite)
            else:
                if not perm_sqlite.exists():
                    sys.exit(f'Unable to find the SQLite database. Please delete or move {self.outdir}, or set "FORCE_OVERWRITE_AND_RE_RUN" to "Yes" in the run.yaml file.')
                else:
                    shutil.copyfile(perm_sqlite, temp_sqlite)
                    pipeline_run(history_file=temp_sqlite)
                    shutil.copyfile(temp_sqlite, perm_sqlite)

            # Print out the pipeline graph
            pipeprintgraph(make_path(self.outdir, "pipeline_graph.svg"), "svg")


if __name__ == "__main__":
    import doctest

    doctest.testmod()
