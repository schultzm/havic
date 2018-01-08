class Check_dependency:
    """
    Check dependency for analysis pipeline.
    """

    def __init__(self, software):
        self.software = software

    def check_software(self):
        """
        Check if software is installed.
        """
        import shutil
        import os
        import sys
        #os.X_OK checks if the file is executable
        path = shutil.which(self.software, mode=os.X_OK)
        if path is not None:
            print(f'{self.software.ljust(28)}: ok ({path})')
        else:
            print(
                f'Dependency {self.software} not callable in path',
                file=sys.stderr)
