#! /usr/bin/env python

#execution script for Snakemake workflow in python

import argparse
import os.path
import snakemake
import sys
import pprint
import json


thisdir = os.path.abspath(os.path.dirname(__file__))


def main(args):
    #find the Snakefile
    snakefile = os.path.join(thisdir, 'Snakefile')
    if not os.path.exists(snakefile):
        sys.stderr.write('Error: cannot find Snakefile at {}\n'.format(snakefile))
        sys.exit(-1)
    
    #run
    status = snakemake.snakemake(snakefile,configfiles=[args.configfile],
                                 targets=[args.workflow],
                                 forceall=args.forceall,
                                 dryrun=args.dry_run, use_conda=True, unlock=args.unlock, printdag=args.dag,  cores=args.cores)

    if status: 
       return 0
    return 1

    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Python executable for Snakemake workflows', usage='''./run.py -w chosen_workflow -c config.json -other_argument -q X
    
     ''')

    parser.add_argument('-n', '--dry-run', action='store_true', help='run Snakemake in dry-run mode')
    parser.add_argument('-u', '--unlock', action='store_true', help='unlock the current working directory')
    parser.add_argument('-d', '--dag', action='store_true', help='print the DAG of jobs in dot graphviz language. Remember to specify "|dot -Tsvg > dag.svg"')
    parser.add_argument('-q', '--cores', action='store', type=int,  help='specify how many cores to use')
    parser.add_argument('-f', '--forceall', action='store_true', help='force the execution of the workflow even if it is done already')
    parser.add_argument('-w', '--workflow', action='store', help='choose the workflow between: "mapping", "calling", "annotating" and "complete"')
    parser.add_argument('-c', '--configfile', action='store', help='choose the configuration file')
    args = parser.parse_args()

    sys.exit(main(args))
