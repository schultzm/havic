#!/usr/bin/env python3

'''
Parse fasta files (remove spaces in descriptor etc)
Concatenate to single file
'''


class Fasta_file:
    def __init__(self, filename):
        '''
        Initialise the class instance with a filename.
        '''
        import os
        if os.path.exists(os.path.abspath(self.filename)):
            self.filename = filename
        else:
            IOerror('File not found')


def main():
    '''
    The main routine.
    '''
    import argparse
    import os
    import sys

    
if __name__ == '__main__':
    main()

