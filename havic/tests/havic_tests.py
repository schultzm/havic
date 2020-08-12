"""
Unit Tests.
    Copyright (C) 2020 Dr Mark B Schultz dr.mark.schultz@gmail.com
    https://github.com/schultzm/HAVIC.git GNU Affero General Public License
    <https://www.gnu.org/licenses/>.
"""


import unittest
import pkg_resources
from random import choice as rndm
from string import ascii_letters
from pathlib import Path, PurePath
from Bio import SeqIO
from .. import (__parent_dir__,
                __ref_seq__,
                __ref_amplicon__,
                __test_seqs__,
                __test_seqs_totrim__,
                __version__
                )
from ..utils.pipeline_runner import Pipeline

PREFIX = f"_{''.join([rndm(ascii_letters) for i in range(4)])}_"
OUTDIR = f"tmpHAVIC_{''.join([rndm(ascii_letters) for i in range(10)])}"

class MergeTestCasePass(unittest.TestCase):
    def setUp(self):
        self.refseq    = SeqIO.read(open(pkg_resources. \
                                    resource_filename(__parent_dir__,
                                                      __ref_seq__), "r"),
                                    "fasta")
        self.testseqs  = list(SeqIO. \
                              parse(open(pkg_resources. \
                                         resource_filename(__parent_dir__,
                                                           __test_seqs__),
                                         "r"), "fasta"))
        self.refamplicon = SeqIO.read(open(pkg_resources. \
                                           resource_filename(__parent_dir__,
                                                             __ref_amplicon__),
                                           "r"), "fasta")
        self.version   = __version__

    def versioner(self):
        """
        Test HAVIC version is not None
        """
        # from .. import __version__
        self.assertFalse(self.version == None)

    def refseqer(self):
        """
        Check refseq id from seq header.
        """
        self.assertEqual(self.refseq.id, 'NC_001489.1')

    def exampler(self):
        """
        Take a peek in the example.fa and see if it is parsing correctly.
        """
        self.assertEqual(self.testseqs[1].id, 'AY644337_55443_seq_1')

    def havnetampliconer(self):
        """
        Parse the reference amplicon, which excludes the primer sites.
        """
        self.assertEqual(self.refamplicon.id, "NC_001489_1_ampliconseq_IB")

    def suite_runner(self):
        """
        Run the pipeline using the full pipeline demo suite.
        """
        detection_pipeline = Pipeline([pkg_resources. \
                                      resource_filename(__parent_dir__,
                                                        __test_seqs__)], #inseqs
                                      __test_seqs_totrim__, #trim to amplicon
                                      None, # if None, use inbuilt refgenome
                                      False, # redo IQTree? No.
                                      3, # How many SNPs?
                                      300, # SNPs in what seq len?
                                      True, #draw the plots
                                      PREFIX, # prepend this to outfiles
                                      OUTDIR, # send files here
                                      5, # k-mer size for minimap2
                                      "ClusterPicker", # tell me where CP is
                                      4) # how many threads to use for IQTree
        for key, value in detection_pipeline.__dict__.items():
            print(f"{key}: {value}\n")
        detection_pipeline._run()


    def pdfs_checker(self):
        """
        Check for two PDF files in OUTDIR
        """
        self.assertTrue(len(list(Path(OUTDIR).glob("*.pdf"))) == 2)

    def csvs_checker(self):
        """
        Check for two CSV files in OUTDIR
        """
        self.assertTrue(len(list(Path(OUTDIR).glob("*.csv"))) == 2)