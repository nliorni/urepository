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


    #find the workflow file

    workflowfile = None
    if os.path.exists(args.workflowfile) and not os.path.isdir(args.workflowfile):
        workflowfile = args.workflowfile
    else:
        for suffix in ('', '.json'):
            tryfile = os.path.join(thisdir, args.workflowfile + suffix)
            if os.path.exists(tryfile) and not os.path.isdir(tryfile):
                sys.stderr.write('Found workflowfile at {}\n'.format(tryfile))
                workflowfile = tryfile
                break

    if not workflowfile:
        sys.stderr.write('Error: cannot find workflowfile {}\n'.format(args.workflowfile))
        sys.exit(-1)

    with open(workflowfile, 'rt') as fp:
        workflow_info = json.load(fp)

    target = workflow_info['workflow_target']

    #l'idea è quello di creare tanti files json quanti tipi di workflow si vogliono utilizzare a partire da una base comune.
    #ad esempio, uno per campione singolo, uno per trio sequencing, uno per rna, ma in cui tutte le regole sono inserite 
    #in un unico workflow

     # run!!
    
    status = snakemake.snakemake(snakefile,#, configfile=paramsfile,
                                 targets=[target],
                                 dryrun=args.dry_run, use_conda=args.use_conda, unlock=args.unlock, printdag=args.dag) #, cores=args.cores)#, config=config)

    if status: # translate "success" into shell exit code of 0
       return 0
    return 1

    #DA FARE: far funzionare l'opzione "cores"
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Python executable for Snakemake workflows', usage='''./run.py <workflow_file> -argument
             ____
            / . .\\                 
            \    ---<
             \  /
   __________/ /
-=:___________/
Snakemake workflows
     ''')

    parser.add_argument('workflowfile')
    parser.add_argument('-n', '--dry-run', action='store_true', help='run Snakemake in dry-run mode')
    parser.add_argument('-c', '--use-conda', action='store_true', help='use Conda environment for wrappers')
    parser.add_argument('-u', '--unlock', action='store_true', help='unlock the current working directory')
    parser.add_argument('-d', '--dag', action='store_true', help='print the DAG of jobs in dot graphviz language. Remember to specify "|dot -Tsvg > dag.svg"')
    #parser.add_argument('-q', '--cores', action='store', help='specify the cores to use')
    args = parser.parse_args()

    sys.exit(main(args))
