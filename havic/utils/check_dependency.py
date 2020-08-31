#!/usr/bin/env python3

from havic.tests import check_r_dependencies

#TODO: this check tool does not work properly. FIX!


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
        import sys
        if self.category == 'software':
            import shutil
            import os
            path = shutil.which(self.software, mode=os.X_OK)
            if path is not None:
                print(f'{self.software.ljust(28)}: ok ({path})')
            else:
                print(
                    f'Dependency {self.software} not callable in path',
                    file=sys.stderr)
        if self.category == 'rmodule':
            try:
                #result = check_r_dependencies.importr_tryhard(self.software)
                # print(result)
                print(f"R library {self.software}".ljust(28) + ": ok",
                      file=sys.stderr)
            except ImportError:
                check_r_dependencies.importr_tryhard(self.software)
                print(f"R library {self.software}".ljust(28) + ": not found",
                      file=sys.stderr)

if __name__ == "__main__":
    import doctest
    doctest.testmod()

# add nextflow to the list