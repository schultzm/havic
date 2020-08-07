"""
Unit Tests.
    Copyright (C) 2020 Dr Mark B Schultz dr.mark.schultz@gmail.com
    https://github.com/schultzm/HAVIC.git GNU Affero General Public License
    <https://www.gnu.org/licenses/>.
"""

import unittest
from .. import (__parent_dir__,
                __ref_seq__,
                __test_seqs__,
                __test_seqs_totrim__,
                __version__
                )
import pkg_resources
from ..utils import *


class MergeTestCasePass(unittest.TestCase):
    def setUp(self):
        self.refseq    = pkg_resources.resource_filename(__parent_dir__,
                                                         __ref_seq__)
        self.version   = __version__

    def refseqer(self):
        """
        Check refseq id from seq header.
        """
        from Bio import SeqIO
        seqobj = SeqIO.read(open(self.refseq, 'r'), 'fasta')
        self.assertEqual(seqobj.id, 'NC_001489.1')

    def versioner(self):
        """
        Test HAVIC version is not None
        """
        # from .. import __version__
        self.assertFalse(self.version == None)

