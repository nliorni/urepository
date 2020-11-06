#! /usr/bin/env python
"""
Execution script for snakemake workflows.
"""
import argparse
import os.path
import snakemake
import sys
import pprint
import json


thisdir = os.path.abspath(os.path.dirname(__file__))


def main(args):
    # first, find the Snakefile
    snakefile = os.path.join(thisdir, 'Snakefile')
    if not os.path.exists(snakefile):
        sys.stderr.write('Error: cannot find Snakefile at {}\n'.format(snakefile))
        sys.exit(-1)

    print('--------')
    print('details!')
    print('\tsnakefile: {}'.format(snakefile))
    print('--------')

     # run!!
    
    status = snakemake.snakemake(snakefile,#, configfile=paramsfile,
                                 #targets=[target], printshellcmds=True,
                                 dryrun=args.dry_run)#, config=config)

    if status: # translate "success" into shell exit code of 0
       return 0
    return 1
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='run snakemake workflows', usage='''run <workflow> <parameters> [<target>]
Run snakemake workflows, using the given workflow name & parameters file.
''')

    #parser.add_argument('workflowfile')
    #parser.add_argument('paramsfile')
    parser.add_argument('-n', '--dry-run', action='store_true', help='run Snakemake in dry-run mode')
    args = parser.parse_args()

    sys.exit(main(args))
