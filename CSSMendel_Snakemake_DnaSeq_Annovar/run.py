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

     # run!!
    status = snakemake.snakemake(snakefile, configfiles=[args.configfile],
                                 targets=[args.workflow],
                                 forceall=args.forceall,
                                 dryrun=args.dry_run, use_conda=True, force_incomplete=args.reruninc, unlock=args.unlock, printdag=args.dag,  cores=args.cores, lint=args.lint)

    if status: 
       return 0
    return 1

    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Python executable for Snakemake workflows', usage='''./run.py -w workflow -other_argument -c config.json(yaml) -q X
             ____
            / . .\\                 
            \    ---< "Snake crawls... Snakemake workflows"
             \  /
   __________/ /
-=:___________/
     ''')

    parser.add_argument('-n', '--dry-run', action='store_true', help='-run Snakemake in dry-run mode')
    parser.add_argument('-ri', '--reruninc', action='store_true', help='-rerun rules with incomplete outputs')
    parser.add_argument('-u', '--unlock', action='store_true', help='-unlock the current working directory')
    parser.add_argument('-d', '--dag', action='store_true', help='-print the DAG of jobs in dot graphviz language. Remember to specify "|dot -Tsvg > dag.svg"')
    parser.add_argument('-q', '--cores', action='store', type=int,  help='-[required] specify how many cores to use')
    parser.add_argument('-f', '--forceall', action='store_true', help='-force the execution of the workflow even if it is done already')
    parser.add_argument('-w', '--workflow', action='store', help='-[required] choose the workflow. Options: gatk, deepvariant, single_gatk, single_deepvariant')
    parser.add_argument('-c', '--configfile', action='store', help='-[required] choose the configuration file')
    parser.add_argument('-l', '--lint', action='store', type=str, help='print the lint. Default none, can be "plenty" or "json"' )
    args = parser.parse_args()

    sys.exit(main(args))
