#!/usr/bin/env python3

from rpy2.robjects.packages import importr, isinstalled
import pandas as pd
import shutil
import os


def importr_tryhard(packname):
    """Check for R library

    Args:
        packname (string): The name of the R package

    Returns:
        class 'SignatureTranslatedPackage' : rpack name
    """
    # print(InstalledPackages())
    rpack = None
    if isinstalled(packname):
        rpack = importr(packname, on_conflict="warn")
        # # cmd = f'echo "library(\'{packname}\'); '
        #         f'path.package(\'{packname}\')" | R --no-save'
    return rpack


class Dependency:
    """
    Check dependency for analysis pipeline.

    >>> Dependency('grep', 'software').check()
    grep                        : ok (/usr/bin/grep)
    """

    def __init__(self, software, category):
        """
        Initialize the class
        :param software:
        :param category:
        >>> x = Dependency('grep', 'software')
        >>> print(x.software)
        grep
        >>> print(x.category)
        software
        """
        self.software = software
        self.category = category.lower()

    def check(self):
        """
        Check if software is installed.
        :return: result to stdout
        >>> Dependency('grep', 'software').check()
        grep                        : ok (/usr/bin/grep)
        >>> Dependency('Rsamtools', 'rmodule').check()
        """
        path = None
        if self.category == "software":
            result = shutil.which(self.software, mode=os.X_OK)
            if result:
                path = pd.DataFrame(
                    {"type": "software", "status": "ok"}, index=[self.software]
                )
            else:
                path = pd.DataFrame(
                    {"type": "software", "status": "not found"}, index=[self.software]
                )
        elif self.category == "rmodule":
            result = importr_tryhard(self.software)
            if result:
                path = pd.DataFrame(
                    {"type": "Rlib", "status": "ok"}, index=[self.software]
                )
            else:
                path = pd.DataFrame(
                    {"type": "Rlib", "status": "not found"}, index=[self.software]
                )
        else:
            pass
        return path


if __name__ == "__main__":
    import doctest

    doctest.testmod()
