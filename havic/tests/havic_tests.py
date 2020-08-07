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

    def versioner(self):
        # from .. import __version__
        self.assertFalse(self.version == None)

