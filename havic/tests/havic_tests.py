"""
Unit Tests.
    Copyright (C) 2020 Dr Mark B Schultz dr.mark.schultz@gmail.com
    https://github.com/schultzm/HAVIC.git GNU Affero General Public License
    <https://www.gnu.org/licenses/>.
"""


import unittest
# import pkg_resources
from pkg_resources import resource_filename as rf
from random import choice as rndm
from string import ascii_letters
from pathlib import Path, PurePath
import yaml
import sys
from Bio import SeqIO
from .. import (__parent_dir__,
                __havic_yaml__,
                __version__
                )
from ..utils.pipeline_runner import Pipeline

class MergeTestCasePass(unittest.TestCase):
    def setUp(self):
        self.version   = __version__
        self.yaml = yaml.load(open(rf(__parent_dir__,
                                   __havic_yaml__)
                                  ),
                             Loader=yaml.FullLoader)
        # print(self.yaml)

    def versioner(self):
        """
        Test HAVIC version is not None
        """
        self.assertFalse(self.version == None)

    def suite_runner(self):
        """
        Run the pipeline using the full pipeline demo suite.
        """
        detection_pipeline = Pipeline(self.yaml)
        # for key, value in detection_pipeline.__dict__.items():
        #     print(f"{key}: {value}\n")
        self.assertEqual(detection_pipeline.yaml_in['QUERY_FILES'][0],
                         'data/example1.fa')
        detection_pipeline._run()


    def pdfs_checker(self):
        """
        Check for two PDF files in OUTDIR
        """
        self.assertTrue(len(list(Path(self.yaml['OUTDIR']).glob("*.pdf"))) >= 2)

    def csvs_checker(self):
        """
        Check for two CSV files in OUTDIR
        """
        self.assertTrue(len(list(Path(self.yaml['OUTDIR']).glob("*.csv"))) == 2)