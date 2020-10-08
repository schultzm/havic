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
    subparser_args2 = argparse.ArgumentParser(add_help=False)
    subparser_args2.add_argument("test_suite", choices=["hav_amplicon",
                                          "hav_wgs",
                                          "measles_wgs",
                                          "hiv_amplicon"],
                                help="""The test suite to run.""")#,
                                # dest="testsuite")

    subparser_modules = parser.add_subparsers(
        title="Sub-commands help", help="", metavar="", dest="subparser_name"
    )
    subparser_modules.add_parser(
        "detect",
        help="""Detect infection clusters from cDNA or 
        DNA consensus sequences.""",
        description="Start the infection cluster detection pipeline.",
        parents=[subparser_args1],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    subparser_modules.add_parser(
        "version", help="Print version.", description="Print version."
    )
    subparser_modules.add_parser(
        "test",
        parents=[subparser_args2],
        help="Run havic test using pre-packaged example data.",
        description="Run havic test using pre-packaged example data.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    args = parser.parse_args()

    if not args.subparser_name:
        parser.print_help()

    if args.subparser_name == "test":
        import unittest
        from .tests.suite_test import suite, suite2, suite3, suite4

        runner = unittest.TextTestRunner(verbosity=2)
        if args.test_suite == "hav_amplicon":
            runner.run(suite())
        elif args.test_suite == "hav_wgs":
            runner.run(suite2())
        elif args.test_suite == "measles_wgs":
            runner.run(suite3())
        else:
            runner.run(suite4()) # the hiv_amplicon test data

    elif args.subparser_name == "detect":
        import yaml

        yaml_in = yaml.load(open(args.yaml_path, "r"), Loader=yaml.FullLoader)
        from .utils.pipeline_runner import Pipeline

        detection_pipeline = Pipeline(yaml_in)
        detection_pipeline._run()
        get_execution_time(yaml_in["OUTDIR"])

    elif args.subparser_name == "version":
        from .utils.version import Version

        Version()


if __name__ == "__main__":
    main()
