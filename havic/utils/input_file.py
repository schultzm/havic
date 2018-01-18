#!/usr/bin/env python3

class Input_file:
    '''

    Input files for analysis.
    '''

    def __init__(self, filename, file_category):
        """
        Initialise the class instance with a filename.
        """
        import os
        if os.path.isfile(os.path.abspath(filename)):
            self.filename = os.path.abspath(filename)
        else:
            raise IOError(f'{file_category} file {filename} not found.')

if __name__ == "__main__":
    import doctest
    doctest.testmod()