#!/usr/bin/env python3
"""
A pipeline for aligning amplicons and picking transmission clusters.

Steps:
    Parse fasta files (remove spaces in descriptor etc).
    Concatenate to single file.
    Run the pipeline.
    Trim the alignment.

todo: doctest in classes
"""

from .utils.input_file import Input_file
from .mapping.bam2fasta import bam2fasta
from .plotters.treeplot_snpplot import plot_functions
from .data.havnet_amplicon import havnet_ampliconseq
from . import __ref_seq__, __parent_dir__
import pkg_resources
import io


def main():
    """Perform the main routine."""
    import argparse
    import os
    import sys
    parser = argparse.ArgumentParser(
        prog='havtrans',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    subparsers = parser.add_subparsers(
        title="Sub-commands help", help="", metavar="", dest="subparser_name")
    subparser = subparsers.add_parser(
        "run", help="Run the analysis.", description="Run the pipeline.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    subparser.add_argument(
        "-q", "--query_files", help="Query file", nargs="+", required=True)
    subparser.add_argument(
        "-t",
        "--trim_seqs",
        help="""Fasta headers of sequences to be trimmed to match length of
                reference amplicon in the alignment.""",
        nargs="+",
        required=False)
    subparser.add_argument(
        "-s",
        "--subject_file",
        help="""Subject file.
                Default is the complete HAVNET reference genome:
                NC_001489.1 Hepatitis A virus.""",
        default=None,
        required=False)
    subparser.add_argument(
        "-r",
        "--redo",
        help="Redo all  (force redo).",
        default=False,
        action="store_true",
        required=False)
    subparser.add_argument(
        "-n",
        "--n_snps",
        help="""Number of SNPS in distance
                                       fraction (numerator).""",
        default=3,
        type=int,
        required=False)
    subparser.add_argument(
        "-l",
        "--seqlen",
        help="""Sequence length in distance
                                       fraction (denominator).""",
        default=300,
        type=int,
        required=False)
    subparser.add_argument(
        "-p",
        "--prefix",
        help="""Filename prefix.""",
        default="_test_",
        required=False)
    subparser.add_argument(
        "-k",
        "--minimap2_kmer",
        help="""k-mer size for minimap2 step.""",
        default=5,
        type=int,
        choices=[3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27],
        required=False)
    subparsers.add_parser(
        "version", help="Print version.", description="Print version.")
    subparsers.add_parser(
        "check", help="Check dependencies are in path.",
        description="Check dependencies.")

    args = parser.parse_args()
    import sys
    if not args.subparser_name:
        parser.print_help()
    elif args.subparser_name == 'check':
        from .utils.check_dependency import Dependency
        from .tests.dependencies import SOFTWAREZ, R_LIBS
        for software in SOFTWAREZ:
            Dependency(software, 'software').check()
        for rlib in R_LIBS:
            Dependency(rlib, 'rmodule').check()
    elif args.subparser_name == 'version':
        from .utils.version import Version
        Version()
    elif args.subparser_name == 'run':
        # 0 Read in the infiles
        queries = [Input_file(file, "Query").filename for file in
                   args.query_files]
        if args.subject_file is None:
            subject = pkg_resources.resource_filename(__parent_dir__,
                                                      __ref_seq__)
        else:
            subject = Input_file(args.subject_file, "Subject").filename
        print(f"Will map amplicons to {subject}")

        # 1 Compile the fasta files to single file
        from Bio import SeqIO
        quality_controlled_seqs = []
        tmp_fasta = os.path.expanduser(f"~/{args.prefix}all_tmp.fa")
        tmp_bam = os.path.expanduser(f"~/{args.prefix}HAV_all_minimap2.bam")
        fasta_from_bam = os.path.expanduser(
            f"~/{args.prefix}HAV_all_minimap2.stack.fa")
        fasta_from_bam_trimmed = os.path.expanduser(
            f"~/{args.prefix}HAV_all_minimap2" + ".stack.trimmed.fa")
        for query_file in queries:
            print(query_file)
            for record in SeqIO.parse(query_file, "fasta"):
                # 1.01 Fix fasta headers
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
        # 1.01 Append the reference amplicon to the alignment
        quality_controlled_seqs.append(
            SeqIO.read(io.StringIO(havnet_ampliconseq), "fasta"))
        SeqIO.write(quality_controlled_seqs, tmp_fasta, "fasta")
        # todo - 1.1 trim the sequences to remove primers
        # 2 get ref and ref stats
        refseq = SeqIO.read(subject, "fasta")
        reflen = len(refseq.seq)
        header = refseq.id
        # 3 get minimap2 done
        import os
        cmd = f"minimap2 -k {args.minimap2_kmer} -a {subject} {tmp_fasta} " \
              f"| samtools sort > {tmp_bam}"
        print(cmd)
        os.system(cmd)
        cmd = f"samtools index {tmp_bam}"
        os.system(cmd)
        # 3.1 print the unmapped sequences.
        cmd = f"samtools view -f 4 {tmp_bam} | cut -f 1"
        print(f"Unmapped reads at k-mer {args.minimap2_kmer}:")
        os.system(cmd)
        # 4 get the fasta from the bam using bam2fasta
        from rpy2 import robjects
        import warnings
        from rpy2.rinterface import RRuntimeWarning
        warnings.filterwarnings("ignore", category=RRuntimeWarning)
        try:
            cmd = bam2fasta % (tmp_bam, f"{tmp_bam}.bai", header, 1, reflen,
                               fasta_from_bam)
            print(cmd)
            robjects.r(cmd)
        except OSError:
            sys.exit("bam2fasta error")
        # 4.0 Import the alignment
        from Bio import AlignIO
        alignment = AlignIO.read(open(fasta_from_bam, "r"), "fasta")
        # 4.1 Trim the alignment for isolates in arg.trim_seq to match
        # refamplicon.
        aln_trim = Trimmed_alignment(alignment,
                                     SeqIO.read(
                                         io.StringIO(havnet_ampliconseq),
                                         "fasta").id, '-', args.trim_seqs)
        aln_trim._get_refseq_boundary()
        aln_trim.trim_seqs_to_ref()
        # 4.1.1 Depad the alignment.
        aln_trim.depad_alignment()
        AlignIO.write(aln_trim.alignment, fasta_from_bam_trimmed, "fasta")
        # 5 Run iqtree on the extracted bam2fasta
        if args.redo:
            redo = "redo -redo"
        else:
            redo = " TEST"
        cmd = f"iqtree -s {fasta_from_bam_trimmed} -nt AUTO -bb 1000 -m{redo}"
        # TN+I+G4{redo}"
        os.system(cmd)
        # 5.1 Midpoint root the phylogeny using ete3
        treefile = f"{fasta_from_bam_trimmed}.treefile"

        print(f"Reading {treefile}")
        mp_treefile = f"{fasta_from_bam_trimmed}.mp.treefile"
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
        cmd = f"java -jar {CLUSTER_PICKER} {fasta_from_bam_trimmed} " \
              f"{mp_treefile} 70 95 {args.n_snps/args.seqlen} 15 valid"
        print(cmd)
        os.system(cmd)
        cmd = f"cp {fasta_from_bam_trimmed}.mp_clusterPicks.nwk.figTree " \
              f"{fasta_from_bam_trimmed}.div_{args.n_snps}SNPsIn{args.seqlen}" \
              f"bp.mp_clusterPicks.nwk.figTree"
        os.system(cmd)
        # 6 Link tree to alignment and plot it
        # treestring = open(f"{fasta_from_bam_trimmed}.mp.treefile", "r").read()
        with open(f"{fasta_from_bam_trimmed}.Rplot.R", "w") as out_r:
            cmd = plot_functions.replace("<- z", "<- \"" +
                                         fasta_from_bam_trimmed + "\"") \
                .replace("<- a", "<- " + str(args.n_snps)) \
                .replace("<- b", "<- " + str(args.seqlen)) \
                .replace("<- k", "<- " + str(args.minimap2_kmer))
            print(cmd)
            out_r.write(cmd)
        os.system(f"Rscript {fasta_from_bam_trimmed}.Rplot.R")


if __name__ == "__main__":
    main()
