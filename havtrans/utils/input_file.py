class Input_file:
    '''

    Input files for analysis.
    '''

    def __init__(self, filename, file_category):
        """
        Initialise the class instance with a filename.

        >>> fname = Input_file('tests/test_headers.fa', 'fasta')
        >>> print(fname.filename)
        /Users/mschultz/tests/test_headers.fa
        """
        import os
        if os.path.exists(os.path.abspath(filename)):
            self.filename = os.path.abspath(filename)
        else:
            raise IOError(f'{file_category} file {filename} not found.')
