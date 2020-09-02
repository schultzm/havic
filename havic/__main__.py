#!/usr/bin/env python3
"""
A pipeline for aligning amplicons and picking transmission clusters.

Steps:
    See the output svg file after running.

todo: doctest in classes
"""

from datetime import datetime
from pathlib import PurePath

STARTTIME = datetime.now()


def get_execution_time(outdir):
    "Calculate length of run."
    print(f"\nTotal runtime (HRS:MIN:SECS): {str(datetime.now() - STARTTIME)}")
    print(f"Results in {PurePath(outdir)} 😷")


def main():
    """Perform the main routine."""
    import argparse

    parser = argparse.ArgumentParser(
        prog="havic", formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    subparser_args1 = argparse.ArgumentParser(add_help=False)
    subparser_args1.add_argument("yaml_path", help="""Path to yaml config.""")
    subparser_modules = parser.add_subparsers(
        title="Sub-commands help", help="", metavar="", dest="subparser_name"
    )
    subparser_modules.add_parser(
        "detect",
        help="""Detect Hepatitis A Virus infection clusters from HAVNET
                protocol amplicon sequences.""",
        description="Start the infection cluster detection pipeline.",
        parents=[subparser_args1],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    subparser_modules.add_parser(
        "version", help="Print version.", description="Print version."
    )
    subparser_modules.add_parser(
        "depcheck",
        help="Check dependencies are in path.  Requires Rpy2.",
        description="Check dependencies.",
    )
    subparser_modules.add_parser(
        "test",
        help="Run HAVIC test using pre-packaged example data.",
        description="Run HAVIC test using pre-packaged example data.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    args = parser.parse_args()

    if not args.subparser_name:
        parser.print_help()
    elif args.subparser_name == "depcheck":
        from .utils.check_dependency import Dependency
        from .data.dependencies import SOFTWAREZ, R_LIBS

        for software in SOFTWAREZ:
            Dependency(software, "software").check()
        for rlib in R_LIBS:
            Dependency(rlib, "rmodule").check()

    elif args.subparser_name == "test":
        import unittest
        from .tests.suite_test import suite

        runner = unittest.TextTestRunner(verbosity=2)
        runner.run(suite())

    elif args.subparser_name == "detect":
        from .utils.pipeline_runner import Pipeline

        detection_pipeline = Pipeline(
            (args.query_fofn, False),  # test=False
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
            args.iqtree_threads,
        )
        for key, value in detection_pipeline.__dict__.items():
            print(f"{key}: {value}\n")
        detection_pipeline._run()
        get_execution_time(args.outdir)

    elif args.subparser_name == "version":
        from .utils.version import Version

        Version()


if __name__ == "__main__":
    main()
