#!/usr/bin/env python3
"""
A pipeline for aligning amplicons and picking transmission clusters.

Steps:
    See the output svg file after running.

todo: doctest in classes
"""

from datetime import datetime
import random
import string

STARTTIME = datetime.now()

def get_execution_time(outdir):
    import os
    "Calculate length of run."
    print(
        f"\nTotal runtime (HRS:MIN:SECS): {str(datetime.now() - STARTTIME)}")
    print(f"Results in {os.path.abspath(outdir)} 😷")


def main():
    """Perform the main routine."""
    import argparse
    parser = argparse.ArgumentParser(
        prog='havic',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    subparser_args1 = argparse.ArgumentParser(add_help=False)
    subparser_args1.add_argument(
        "-t",
        "--trim_seqs",
        help="""Fasta headers of sequences to be trimmed to match length of
                reference amplicon in the alignment.""",
        nargs="+",
        default='',
        required=False)
    subparser_args1.add_argument(
        "-s",
        "--subject_file",
        help="""Subject file.  If not provided,
                uses the pre-packaged complete HAVNET reference genome:
                NC_001489.1 Hepatitis A virus.""",
        default=None,
        required=False)
    subparser_args1.add_argument(
        "-r",
        "--redo",
        help="Redo all  (not yet implemented).",
        default=False,
        action="store_true",
        required=False)
    subparser_args1.add_argument(
        "-n",
        "--n_snps",
        help="""Number of SNPS in distance
                fraction (numerator) for ClusterPicker assingments.""",
        default=3,
        type=int,
        required=False)
    subparser_args1.add_argument(
        "-l",
        "--seqlen",
        help="""Sequence length in distance
                fraction (denominator) for ClusterPicker assignments.""",
        default=300,
        type=int,
        required=False)
    subparser_args1.add_argument(
        "-m",
        "--matrixplots",
        help="""Switch on PDF plotting of: Multiple Sequence Alignment (MSA)
                next to IQ-Tree Maximum Likelihood tree), and; heatmap of SNP
                distances.""",
        default=False,
        action="store_true",
        required=False)
    subparser_args1.add_argument(
        "-p",
        "--prefix",
        help="""Filename prefix.""",
        default=f"_{''.join([random.choice(string.ascii_letters) for i in range(4)])}_",
        required=False)
    subparser_args1.add_argument(
        "-o",
        "--outdir",
        help="""Output directory.""",
        default=f"tmpHAVIC_{''.join([random.choice(string.ascii_letters) for i in range(10)])}",
        required=False)
    subparser_args1.add_argument(
        "-k",
        "--minimap2_kmer",
        help="""k-mer size for minimap2 step.""",
        default=5,
        type=int,
        choices=[3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27],
        required=False)
    import os
    subparser_args1.add_argument(
        "-c",
        "--path_to_clusterpicker",
        help="""Path to ClusterPicker file.  Follow instructions at
                http://hiv.bio.ed.ac.uk/software.html""",
        default="ClusterPicker",
        required=False)
    from multiprocessing import cpu_count
    subparser_args1.add_argument(
        "-z",
        "--iqtree_threads",
        help="""Number of threads for IQtree""",
        default=4,
        type=int,
        choices=[2, 4],
        required=False)
    subparser_args2 = argparse.ArgumentParser(add_help=False)
    subparser_args2.add_argument(
        "-q", "--query_files", help="Query file", nargs="+", required=True)

    subparser_modules = parser.add_subparsers(
        title="Sub-commands help", help="", metavar="", dest="subparser_name")
    subparser_modules.add_parser(
        "detect",
        help="""Detect Hepatitis A Virus infection clusters from HAVNET
                protocol amplicon sequences.""",
        description="Start the infection cluster detection pipeline.",
        parents=[subparser_args1, subparser_args2],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    subparser_modules.add_parser(
        "version", help="Print version.", description="Print version.")
    subparser_modules.add_parser(
        "depcheck", help="Check dependencies are in path.  Requires Rpy2.",
        description="Check dependencies.")
    subparser_modules.add_parser(
        "test", help="Run HAVIC test using pre-packaged example data.",
        description="Run HAVIC test using pre-packaged example data.",
        parents=[subparser_args1],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    args = parser.parse_args()

    if not args.subparser_name:
        parser.print_help()
    elif args.subparser_name == 'depcheck':
        try:
            import rpy2
        except ImportError as error:
            import sys
            sys.exit(f'\'havic {args.subparser_name}\' requires rpy2.\n'
                     f'{error.__class__.__name__}\n'
                     f'Please install rpy2 to use this feature.')
        from .utils.check_dependency import Dependency
        from .tests.dependencies import SOFTWAREZ, R_LIBS
        for software in SOFTWAREZ:
            Dependency(software, 'software').check()
    elif args.subparser_name == 'test':
        print(f"Running test suite...")
        from . import __parent_dir__, __test_seqs__, __test_seqs_totrim__
        import pkg_resources
        from .utils.pipeline_runner import Pipeline
        test_query = [pkg_resources.resource_filename(__parent_dir__,
                                                      __test_seqs__)]
        test_seqs_totrim = [test_seq_totrim for test_seq_totrim in
                            __test_seqs_totrim__]
        detection_pipeline = Pipeline(test_query,
                                      test_seqs_totrim,
                                      args.subject_file,
                                      args.redo,
                                      args.n_snps,
                                      args.seqlen,
                                      args.matrixplots,
                                      args.prefix,
                                      args.outdir,
                                      args.minimap2_kmer,
                                      args.path_to_clusterpicker,
                                      args.iqtree_threads)
        for key, value in detection_pipeline.__dict__.items():
            print(f"{key}: {value}\n")
        detection_pipeline._run()
        get_execution_time(args.outdir)

    elif args.subparser_name == 'detect':
        from .utils.pipeline_runner import Pipeline
        detection_pipeline = Pipeline(args.query_files,
                                      args.trim_seqs,
                                      args.subject_file,
                                      args.redo,
                                      args.n_snps,
                                      args.seqlen,
                                      args.matrixplots,
                                      args.prefix,
                                      args.outdir,
                                      args.minimap2_kmer,
                                      args.path_to_clusterpicker,
                                      args.iqtree_threads)
        for key, value in detection_pipeline.__dict__.items():
            print(f"{key}: {value}\n")
        detection_pipeline._run()
        get_execution_time(args.outdir)

    elif args.subparser_name == 'version':
        from .utils.version import Version
        Version()


if __name__ == "__main__":
    main()
