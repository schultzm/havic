#!/usr/bin/env python3
"""
A pipeline for aligning amplicons and picking transmission clusters.

Steps:
    See the output svg file after running.

todo: doctest in classes
"""


import tempfile


def main():
    """Perform the main routine."""
    import argparse
    import os
    import sys
    parser = argparse.ArgumentParser(
        prog='havic',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    subparsers = parser.add_subparsers(
        title="Sub-commands help", help="", metavar="", dest="subparser_name")
    subparser = subparsers.add_parser(
        "detect", help="Start the Hepatitis A Virus infection cluster detection.", description="Start the Hepatitis A Virus infection cluster detection.",
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
        "-o",
        "--outdir",
        help="""Output directory.""",
        default=tempfile.mkdtemp(dir=os.path.abspath('.')),
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
    elif args.subparser_name == 'detect':
        from .utils.pipeline_runner import Pipeline
        detection_pipeline = Pipeline(args.query_files,
                            args.trim_seqs,
                            args.subject_file,
                            args.redo,
                            args.n_snps,
                            args.seqlen,
                            args.prefix,
                            args.outdir,
                            args.minimap2_kmer)
        for key, value in detection_pipeline.__dict__.items():
            print(f"{key}: {value}\n")
        detection_pipeline.run()



if __name__ == "__main__":
    main()
