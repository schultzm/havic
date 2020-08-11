"""
Unit Test suite builder.
    Copyright (C) 2020 Dr Mark B Schultz dr.mark.schultz@gmail.com
    https://github.com/schultzm/gnb.git GNU Affero General Public License
    <https://www.gnu.org/licenses/>.
"""

import unittest
from ..tests.havic_tests import MergeTestCasePass


def suite():
    """
    This is the test suite.
    """
    suite = unittest.TestSuite()
    suite.addTest(MergeTestCasePass("versioner"))
    suite.addTest(MergeTestCasePass("refseqer"))
    suite.addTest(MergeTestCasePass("exampler"))
    return suite