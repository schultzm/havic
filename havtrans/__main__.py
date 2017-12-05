#!/usr/bin/env python3

'''
Parse fasta files (remove spaces in descriptor etc)
Concatenate to single file
'''

from havtrans.classes.classes import Input_file, Check_dependency
from havtrans.tests import check_r_dependencies
from havtrans.tests.dependencies import SOFTWAREZ, R_LIBS, CLUSTER_PICKER
import sys

def main():
    '''
    The main routine.
    '''
    import argparse
    import os
    import sys
    parser = argparse.ArgumentParser(description='Run HAVTrans')
    subparsers = parser.add_subparsers(title='Sub-commands help',
                                       help='', metavar='',
                                       dest='subparser_name')
    subparser_run = subparsers.add_parser('run', help='Run the analysis.',
                                           description='Run the pipeline.')
    subparser_run.add_argument('-q', '--query_files', help='Query file',
                               nargs='+', required=True)
    subparser_run.add_argument('-s', '--subject_file', help='Subject file',
                               required=True)
    subparser_version = subparsers.add_parser('version', help='Print version.',
                                           description='Print version.')
    args = parser.parse_args()
    
    if not args.subparser_name:
        os.system('havtrans -h')
    queries = [Input_file(file, 'Query').filename for file in args.query_files]
    subject = Input_file(args.subject_file, 'Subject').filename
    
    for dep in SOFTWAREZ:
        path = Check_dependency(dep)
        path.check_software()
    for dep in R_LIBS: #move this to class
        check_r_dependencies.importr_tryhard(dep)
        print(f'R library {dep}'.ljust(28)+': ok', file=sys.stderr)


    
if __name__ == '__main__':
    main()

