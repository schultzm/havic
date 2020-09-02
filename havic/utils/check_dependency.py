#!/usr/bin/env python3

from rpy2.robjects.packages import importr, isinstalled

# import rpy2.rinterface as rinterface

# from rpy2.rinterface import RRuntimeWarning
# from rpy2.robjects.help.Package


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
        if self.category == "software":
            import shutil
            import os

            path = shutil.which(self.software, mode=os.X_OK)
            if path is not None:
                print(f"{self.software.ljust(40)}: ok")
            else:
                print(f"Dependency {self.software} not found")
        if self.category == "rmodule":
            rlib = importr_tryhard(self.software)
            if rlib:
                print(f"Rlib {rlib.__name__}".ljust(40) + ": ok")
            else:
                print(f"R library {self.software}".ljust(40) + ": not found")


if __name__ == "__main__":
    import doctest

    doctest.testmod()
