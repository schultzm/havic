#!/usr/bin/env python3
"""
A pipeline for aligning amplicons and picking transmission clusters.

Steps:
    Parse fasta files (remove spaces in descriptor etc).
    Concatenate to single file.
    Run the pipeline.
    Trim the alignment.
"""
import sys
from havtrans.utils.input_file import Input_file
from havtrans.utils.check_dependency import Check_dependency
from havtrans.utils.trim_alignment import Trimmed_alignment
from havtrans.tests import check_r_dependencies
from havtrans.tests.dependencies import SOFTWAREZ, R_LIBS, CLUSTER_PICKER
from havtrans.mapping.bam2fasta import bam2fasta
# from havtrans.plotters.plottree import plottree
from havtrans.plotters.treeplot_snpplot import plot_functions
# from havtrans.plotters.pdfloop import looper
from havtrans.data.havnet_amplicon import havnet_ampliconseq
from havtrans import (__ref_seq__, __parent_dir__, __version__,
                      __version_date__, __author__, __author_email__,
                      __github_username__, __download_url__)
import pkg_resources
import io


def main():
    """Perform the main routine."""
    import argparse
    import os
    import sys
    parser = argparse.ArgumentParser(description="Run HAVTrans")
    subparsers = parser.add_subparsers(
        title="Sub-commands help", help="", metavar="", dest="subparser_name")
    subparser_run = subparsers.add_parser(
        "run", help="Run the analysis.", description="Run the pipeline.")
    subparser_run.add_argument(
        "-q", "--query_files", help="Query file", nargs="+", required=True)
    subparser_run.add_argument(
        "-t",
        "--trim_seqs",
        help="""Fasta headers of sequences to be trimmed to match length of 
                reference amplicon in the alignment.""",
        nargs="+", required=False)
    subparser_run.add_argument(
        "-s",
        "--subject_file",
        help="""Subject file.
               Default is the complete HAVNET reference genome:
               NC_001489.1 Hepatitis A virus.""",
        default=None,
        required=False)
    subparser_run.add_argument(
        "-r",
        "--redo",
        help="Redo all  (force redo).",
        default=False,
        action="store_true",
        required=False)
    subparser_run.add_argument(
        "-n",
        "--n_snps",
        help="""Number of SNPS in distance
                                       fraction (numerator, default=3).""",
        default=3,
        type=int,
        required=False)
    subparser_run.add_argument(
        "-l",
        "--seqlen",
        help="""Sequence length in distance
                                       fraction (denominator, default=300).""",
        default=300,
        type=int,
        required=False)
    subparser_run.add_argument(
        "-p",
        "--prefix",
        help="""Filename prefix.""",
        default="_test_",
        required=False)
    subparser_run.add_argument(
        "-k",
        "--minimap2_kmer",
        help="""k-mer size for minimap2 step.
                                       Default=5.""",
        default=5,
        type=int,
        choices=[3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27],
        required=False)
    subparser_version = subparsers.add_parser(
        "version", help="Print version.", description="Print version.")
    args = parser.parse_args()
    if not args.subparser_name:
        os.system("havtrans -h")
        sys.exit()
    queries = [Input_file(file, "Query").filename for file in args.query_files]

    print(args.trim_seqs, queries)
    if args.subject_file is not None:
        subject = Input_file(args.subject_file, "Subject").filename
    else:
        subject = pkg_resources.resource_filename(__parent_dir__, __ref_seq__)
    print(f"Will map amplicons to {subject}")
    # for dep in SOFTWAREZ:
    #     path = Check_dependency(dep)
    #     path.check_software()
    # for dep in R_LIBS:  # move this to class
    #     check_r_dependencies.importr_tryhard(dep)
    #     print(f"R library {dep}".ljust(28) + ": ok", file=sys.stderr)
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
        SeqIO.read(io.StringIO(havnet_ampliconseq), "fasta"))
    SeqIO.write(quality_controlled_seqs, tmp_fasta, "fasta")
    # 1.1 trim the sequences to remove primers - todo
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
    # 3.1 find the unmapped sequences.
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
    # 4.1 Trim the alignment for isolates in arg.trim_seq to match refamplicon.
    aln_trim = Trimmed_alignment(alignment,
                                 SeqIO.read(
                                     io.StringIO(havnet_ampliconseq),
                                     "fasta").id,
                                 '-', args.trim_seqs)  #._get_refseq_boundary()
    aln_trim._get_refseq_boundary()
    aln_trim.trim_seqs_to_ref()
    # 4.1.1 Convert the Trimmed_alignment object back to instance of MSA.
    aln_trim.depad_alignment()
    print(aln_trim.alignment)
    sys.exit()
    # from Bio.Align import MultipleSeqAlignment
    # alignment = MultipleSeqAlignment(aln_trim.alignment)

    # # 4.2 Trim the whole alignment to get rid of gap-only
    # # sites at 5" and 3" end of aln
    # site_set = {"-"}
    # start_pos = 0
    # while len(site_set) == 1:
    #     site_set = set(alignment[:, start_pos])
    #     start_pos += 1
    # site_set = {"-"}
    # # subtract one due to 0 and 1-based indexing issues
    # end_pos = alignment.get_alignment_length() - 1
    # while len(site_set) == 1:
    #     site_set = set(alignment[:, end_pos])
    #     end_pos -= 1
    # # .format("fasta") #Add 2, again due to indexing discrepancies
    # alignment_trimmed = alignment[:, start_pos:end_pos + 2]
    AlignIO.write(aln_trim.depad_alignment, fasta_from_bam_trimmed, "fasta")
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
