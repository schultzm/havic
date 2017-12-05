class Input_file:
    '''
    Input files for analysis.
    '''
    def __init__(self, filename, file_category):
        '''
        Initialise the class instance with a filename.
        >>> fname = Fasta_file('tests/test_headers.fa')
        >>> print(fname.filename)
        tests/test_headers.fa
        '''
        import os
        if os.path.exists(os.path.abspath(filename)):
            self.filename = os.path.abspath(filename)
        else:
            raise IOError(f'{file_category} file {filename} not found.')


class Check_dependency:
    '''
    Check dependency for analysis pipeline.
    '''
    def __init__(self, software):
        self.software = software
    def check_software(self):
        '''
        Check if software is installed.
        '''
        import shutil
        import os
        import sys
        #os.X_OK checks if the file is executable
        path = shutil.which(self.software, mode=os.X_OK)
        if path is not None:
            print(f'{self.software.ljust(28)}: ok ({path})')
        else:
            print(f'Dependency {self.software} not callable in path',
                  file=sys.stderr)
