"""
Unit Tests.
    Copyright (C) 2020 Dr Mark B Schultz dr.mark.schultz@gmail.com
    https://github.com/schultzm/HAVIC.git GNU Affero General Public License
    <https://www.gnu.org/licenses/>.
"""

import unittest
from .. import (__parent_dir__,
                __ref_seq__,
                __ref_amplicon__,
                __test_seqs__,
                __test_seqs_totrim__,
                __version__
                )
import pkg_resources
from ..utils import *
from Bio import SeqIO


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


        # print(f"Running test suite...")
        # from . import __parent_dir__, __test_seqs__, __test_seqs_totrim__
        # import pkg_resources
        # from .utils.pipeline_runner import Pipeline
        # test_query = [pkg_resources.resource_filename(__parent_dir__,
        #                                               __test_seqs__)]
        # test_seqs_totrim = [test_seq_totrim for test_seq_totrim in
        #                     __test_seqs_totrim__]
        # detection_pipeline = Pipeline(test_query,
        #                               test_seqs_totrim,
        #                               args.subject_file,
        #                               args.redo,
        #                               args.n_snps,
        #                               args.seqlen,
        #                               args.matrixplots,
        #                               args.prefix,
        #                               args.outdir,
        #                               args.minimap2_kmer,
        #                               args.path_to_clusterpicker,
        #                               args.iqtree_threads)
        # for key, value in detection_pipeline.__dict__.items():
        #     print(f"{key}: {value}\n")
        # detection_pipeline._run()
        # get_execution_time(args.outdir)
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
        Parse the string object as a Seq object.  Use the string 
        as opposed to fasta so end user can read the comment that primers are
        excluded from the amplicon.
        """
        self.assertEqual(self.refamplicon.id, "NC_001489_1_ampliconseq_IB")

