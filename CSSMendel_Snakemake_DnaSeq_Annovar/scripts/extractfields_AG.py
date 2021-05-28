#! /usr/bin/env python

import os.path
import os
import argparse


#cerca il fieldsfile

thisdir = os.path.abspath(os.path.dirname(__file__))

fieldsfile = os.path.join(thisdir, 'fieldsfile.txt')
if not os.path.exists(fieldsfile):
    sys.stderr.write('Error: cannot find the fieldsfile in the current folder')
    sys.exit(-1)

    
#introduco command line argument per l'input vcf file
    
parser = argparse.ArgumentParser(description='Specify the input file to extract fields', usage='python3 extractfields.py --input input.vcf --path snpsiftpath > output.csv')
parser.add_argument('-p', '--path', action='store', type=str, help='-specify your SnpSift ExtractFields Path (/path/to/SnpSift extractFields)')
parser.add_argument('-i', '--input', action='store', type=str,  help='-[required] specify the input file')
args = parser.parse_args()

inputfile=args.input
snpsiftef=args.path
#inizializza vettore dei campi

fields_list=[]

#legge il file e inserisce ogni campo nel vettore

with open(fieldsfile, "r") as f:
    for line in f:
       fields_list.append(line.rstrip())


#inizializza il comando da eseguire
#GENERALIZZARE IL PATH DI SNPSIFT

shellcommand= snpsiftef+ ' extractFields -s "," -e "." -info '+inputfile+' '

#completa il comando con i campi di interesse

for i in range(len(fields_list)):
    shellcommand=shellcommand+fields_list[i]+" "

#esegue il comando

os.system(shellcommand)