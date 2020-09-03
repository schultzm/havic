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
        "test",
        help="Run HAVIC test using pre-packaged example data.",
        description="Run HAVIC test using pre-packaged example data.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    args = parser.parse_args()

    if not args.subparser_name:
        parser.print_help()

    if args.subparser_name == "test":
        import unittest
        from .tests.suite_test import suite

        runner = unittest.TextTestRunner(verbosity=2)
        runner.run(suite())

    elif args.subparser_name == "detect":
        import yaml

        yaml_in = yaml.load(open(args.yaml_path, "r"), Loader=yaml.FullLoader)
        from .utils.pipeline_runner import Pipeline

        detection_pipeline = Pipeline(yaml_in)
        # for key, value in detection_pipeline.__dict__.items():
        #     print(f"{key}: {value}\n")
        detection_pipeline._run()
        get_execution_time(yaml_in["OUTDIR"])

    elif args.subparser_name == "version":
        from .utils.version import Version

        Version()


if __name__ == "__main__":
    main()
