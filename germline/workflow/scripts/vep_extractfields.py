#! /usr/bin/env python

import os.path
import os
import argparse
import sys

#introduco command line arguments

parser = argparse.ArgumentParser(description='Specify the input file to extract fields', usage='python3 extractfields.py --input input.vcf --tsv input.tsv')
parser.add_argument('-i', '--input', action='store', type=str,  help='-[required] specify the input file')
parser.add_argument('-f', '--fieldsfile', action='store', type=str, help='-[required] specify the file with the fields to extract')
parser.add_argument('-t', '--tsv', action='store', type=str,  help='-[required] specify the input tsv file to which add annotations')
args = parser.parse_args()

inputfile=args.input
tsvfile=args.tsv
fieldsfile=args.fieldsfile

#inizializza vettore dei campi

fields_list=[]

#legge il file e inserisce ogni campo nel vettore

with open(fieldsfile, "r") as f:
    for line in f:
       fields_list.append(line.rstrip())

#inizializza il comando da eseguire

shellcommand = 'python3 workflow/scripts/vep_annotation_reporter.py -t '+tsvfile+' '+inputfile+' '

#completa il comando con i campi di interesse

for i in range(len(fields_list)):
    shellcommand=shellcommand+fields_list[i]+" "

#esegue il comando

os.system(shellcommand)