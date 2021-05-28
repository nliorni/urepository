import pandas as pd
import argparse

#command line option che chiede il csv in entrata

parser = argparse.ArgumentParser(description='Specify the input csvfile and the output excel file', usage='python3 genelists.py --input file.csv --output file.xlsx')
parser.add_argument('-i', '--input', action='store', type=str, help='-[required] specify the input csvfile')
parser.add_argument('-o', '--output', action='store', type=str, help='[required] specify the output excel file name')

csvfile=args.input
excelfile=args.output

df = pd.read_csv(csvfile, sep='\t')
df.to_excel(excelfile)