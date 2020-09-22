"""
Unit Tests.
    Copyright (C) 2020 Dr Mark B Schultz dr.mark.schultz@gmail.com
    https://github.com/schultzm/HAVIC.git GNU Affero General Public License
    <https://www.gnu.org/licenses/>.
"""


import unittest

from pkg_resources import resource_filename as rf
from pathlib import Path
import yaml
import pandas as pd
import sys
from .. import (__parent_dir__,
                __havic_yaml__,
                __havic_wgs_yaml__,
                __measles_wgs_yaml__,
                __version__)
from ..utils.pipeline_runner import Pipeline
from ..utils.check_dependency import Dependency
from ..data.dependencies import SOFTWAREZ, R_LIBS


class MergeTestCasePass(unittest.TestCase):
    def setUp(self):
        self.version = __version__
        self.yaml = yaml.load(
            open(rf(__parent_dir__, __havic_yaml__)), Loader=yaml.FullLoader
        )
        self.wgsyaml = yaml.load(
            open(rf(__parent_dir__, __havic_wgs_yaml__)), Loader=yaml.FullLoader
        )
        self.measlesyaml = yaml.load(open(rf(__parent_dir__,
                                             __measles_wgs_yaml__)),
                                             Loader=yaml.FullLoader
        )

    def yamler(self):
        """Check yaml loader."""
        self.assertEqual(self.yaml["SUBJECT_FILE"], "data/NC_001489.fa")

    def dependency_checker(self):
        """Check software dependencies."""
        df_list = []
        for software in SOFTWAREZ:
            df_list.append(Dependency(software, "software").check())
        for rlib in R_LIBS:
            df_list.append(Dependency(rlib, "rmodule").check())
        df = pd.concat(df_list)
        df.columns = ["type", "status"]
        print("\n", df, "\n")
        status_set = df.status.unique().tolist()
        if "not found" in status_set:
            sys.exit("Fix software dependencies.  Exiting.")
        self.assertFalse("not found" in status_set)

    def versioner(self):
        """
        Test HAVIC version is not None.
        """
        self.assertTrue(self.version is not None)

    def suite_runner(self):
        """
        Run the pipeline using the full pipeline demo suite.
        """
        detection_pipeline = Pipeline(self.yaml)
        # for key, value in detection_pipeline.__dict__.items():
        #     print(f"{key}: {value}\n")
        self.assertEqual(
            detection_pipeline.yaml_in["QUERY_FILES"][0], "data/example1.fa"
        )
        detection_pipeline._run()

    def pdfs_checker(self):
        """
        Check for two PDF files in OUTDIR.
        """
        pathlist = list(Path(self.yaml["OUTDIR"]).glob("*.pdf"))
        self.assertTrue(len(pathlist) >= 2)

    def csvs_checker(self):
        """
        Check for two CSV files in OUTDIR.
        """
        pathlist = list(Path(self.yaml["OUTDIR"]).glob("*.csv"))
        self.assertTrue(len(pathlist) == 2)

    def svg_checker(self):
        """
        Graphviz functionality gives svg file.
        """
        pathlist = list(Path(self.yaml["OUTDIR"]).glob("*.svg"))
        self.assertTrue(len(pathlist) == 1)

    def wgs_suite_runner(self):
        """
        Run the pipeline using the full pipeline demo suite.
        """
        detection_pipeline = Pipeline(self.wgsyaml)
        # for key, value in detection_pipeline.__dict__.items():
        #     print(f"{key}: {value}\n")
        self.assertEqual(
            detection_pipeline.yaml_in["QUERY_FILES"][0], "data/wgs_hav_seqs/KP879216.fa"
        )
        detection_pipeline._run()

    def measles_suite_runner(self):
        """
        Run the pipeline using the full pipeline demo suite.
        """
        detection_pipeline = Pipeline(self.measlesyaml)
        # for key, value in detection_pipeline.__dict__.items():
        #     print(f"{key}: {value}\n")
        self.assertEqual(
            detection_pipeline.yaml_in["QUERY_FILES"][0], "data/wgs_measles_seqs/measles.fasta"
        )
        detection_pipeline._run()
