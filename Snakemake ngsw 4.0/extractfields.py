#! /usr/bin/env python

import os.path
import os

#cerca il fieldsfile

thisdir = os.path.abspath(os.path.dirname(__file__))

fieldsfile = os.path.join(thisdir, 'fieldsfile.txt')
if not os.path.exists(fieldsfile):
    sys.stderr.write('Error: cannot find the fieldsfile in the current folder')
    sys.exit(-1)

#inizializza vettore dei campi

fields_list=[]

#legge il file e inserisce ogni campo nel vettore

with open(fieldsfile, "r") as f:
    for line in f:
       print(line) 
       fields_list.append(line.rstrip())

print(fields_list)

#inizializza il comando da eseguire

shellcommand='/home/rasrafmekerk/anaconda3/envs/Snakemake/bin/SnpSift extractFields -s "," -e "." -info all.vcf '
shellcommand2="ls -l" #comando di prova

#completa il comando con i campi di interesse

for i in range(len(fields_list)):
    shellcommand=shellcommand+fields_list[i]+" "

print(shellcommand) #vediamo se la stringa è corretta

#esegue il comando

os.system(shellcommand) #esegui comando di prova