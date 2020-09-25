"""
Unit Test suite builder.
    Copyright (C) 2020 Dr Mark B Schultz dr.mark.schultz@gmail.com
    https://github.com/schultzm/gnb.git GNU Affero General Public License
    <https://www.gnu.org/licenses/>.
"""

import unittest
from ..tests.havic_test import (HavAmpliconTestCase,
                               HavWgsTestCase, 
                               MeaslesAmpliconTestCase)


def suite():
    """
    This is the hav amplicon test suite.
    """
    suite_ = unittest.TestSuite()
    suite_.addTest(HavAmpliconTestCase("versioner"))
    suite_.addTest(HavAmpliconTestCase("yamler"))
    # suite_.addTest(HavAmpliconTestCase("dependency_checker"))
    suite_.addTest(HavAmpliconTestCase("suite_runner"))
    suite_.addTest(HavAmpliconTestCase("csvs_checker"))
    suite_.addTest(HavAmpliconTestCase("svg_checker"))

    return suite_

def suite2():
    """
    This is the hav wgs test suite.
    """
    suite_ = unittest.TestSuite()
    suite_.addTest(HavWgsTestCase("wgs_suite_runner"))
    return suite_

def suite3():
    """
    This is the measles amplicon test suite.
    """
    suite_ = unittest.TestSuite()
    suite_.addTest(MeaslesAmpliconTestCase("measles_suite_runner"))
    return suite_