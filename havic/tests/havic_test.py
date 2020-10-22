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
# import pandas as pd
import sys
from .. import (__parent_dir__,
                __havic_yaml__,
                __havic_wgs_yaml__,
                __havic_PMC7259881__,
                __measles_wgs_yaml__,
                __hiv_amplicon_yaml__,
                __version__)
from ..utils.pipeline_runner import Pipeline


class HavAmpliconTestCase(unittest.TestCase):
    def setUp(self):
        self.version = __version__
        self.yaml = yaml.load(open(rf(__parent_dir__, __havic_yaml__)), Loader=yaml.FullLoader)

    def yamler(self):
        """Check yaml loader."""
        self.assertEqual(self.yaml["SUBJECT_FILE"], "data/NC_001489.fa")

    def versioner(self):
        """
        Test HAVIC version is not None.
        """
        self.assertTrue(self.version is not None)

    def suite_runner(self):
        """
        Run the pipeline using the HAV amplicon demo suite.
        """
        detection_pipeline_hav_amplicon = Pipeline(self.yaml)
        detection_pipeline_hav_amplicon._run()
        self.assertTrue(len(list(Path(detection_pipeline_hav_amplicon.outdir).glob("*.pdf"))) >= 2)

    def csvs_checker(self):
        """
        Check for the presence of csv files.
        """
        self.assertTrue(len(list(Path(self.yaml["OUTDIR"]).glob("*.csv"))) >= 2)

    def svg_checker(self):
        """
        Check for the presence of svg file.
        """
        self.assertTrue(len(list(Path(self.yaml["OUTDIR"]).glob("*.svg"))) >= 1)

class HavWgsTestCase(unittest.TestCase):
    def setUp(self):
        self.wgsyaml = yaml.load(open(rf(__parent_dir__, __havic_wgs_yaml__)), Loader=yaml.FullLoader)
        self.detection_pipeline_wgs = Pipeline(self.wgsyaml)

    def wgs_suite_runner(self):
        """
        Run the pipeline using the HAV WGS demo suite.
        """
        self.detection_pipeline_wgs._run()
        self.assertTrue(len(list(Path(self.detection_pipeline_wgs.outdir).glob("*.pdf"))) >= 2)

class MeaslesAmpliconTestCase(unittest.TestCase):
    def setUp(self):
        self.measlesyaml = yaml.load(open(rf(__parent_dir__, __measles_wgs_yaml__)), Loader=yaml.FullLoader)
        self.detection_pipeline_measles = Pipeline(self.measlesyaml)

    def measles_suite_runner(self):
        """
        Run the pipeline using the measles WGS demo suite.
        """
        self.detection_pipeline_measles._run()
        self.assertTrue(len(list(Path(self.detection_pipeline_measles.outdir).glob("*.pdf"))) >= 2)

class HivAmpliconTestCase(unittest.TestCase):
    def setUp(self):
        self.hivyaml = yaml.load(open(rf(__parent_dir__, __hiv_amplicon_yaml__)), Loader=yaml.FullLoader)
        self.detection_pipeline_hiv = Pipeline(self.hivyaml)

    def hiv_suite_runner(self):
        """
        Run the pipeline using the hiv amplicon ClusterPicker demo data.
        """
        self.detection_pipeline_hiv._run()
        self.assertTrue(len(list(Path(self.detection_pipeline_hiv.outdir).glob("*.pdf"))) >= 2)

class HavPmcTestCase(unittest.TestCase):
    def setUp(self):
        self.PMC7259881yaml = yaml.load(open(rf(__parent_dir__, __havic_PMC7259881__)), Loader=yaml.FullLoader)
        self.detection_pipeline_PMC7259881 = Pipeline(self.PMC7259881yaml)

    def pmc_suite_runner(self):
        """
        Run the pipeline using the data from publication PMC7259881.
        """
        self.detection_pipeline_PMC7259881._run()
        self.assertTrue(len(list(Path(self.detection_pipeline_PMC7259881.outdir).glob("*.pdf"))) >= 2)
