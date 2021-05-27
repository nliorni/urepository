import os.path
import argparse
import pandas as pd

thisdir = os.path.abspath(os.path.dirname(__file__))

#inserire command line argument per l'input csvfile:

parser = argparse.ArgumentParser(description='Specify the input csvfile and the gene lists', usage='python3 genelists.py --input file.csv --genelists gl1.txt gl2.txt')
parser.add_argument('-i', '--input', action='store', type=str, help='-[required] specify the input csvfile')
parser.add_argument('-g', '--genelists', action='store', type=str, nargs='+' , help='-[at least one required] specify the gene lists you want to check')
args = parser.parse_args()

csvfile=args.input
genelist=args.genelists

#apri il file csv
#leggine il contenuto

df = pd.read_csv(csvfile, sep='\t')

#per ogni gene list crea una colonna nel file csv di input il cui nome è uguale a quello della gene list

for i in range(len(genelist)):
   listname=genelist[i]
   df[listname]="-"
   df.to_csv(csvfile, sep='\t')

genes=df["ANN[*].GENE"]

#esegui l'attività vera e propria:

for j in range(len(genelist)):
    listname=genelist[j]
    with open(listname, "r") as f:
        for line in f:
            for i in range(len(genes)):
                genarray=genes[i].split(',')
                if line.rstrip('\n') in genarray:
                    df.loc[i,listname]='X'
df.to_csv(csvfile, sep='\t')




