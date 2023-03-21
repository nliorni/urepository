import pandas as pd
import argparse

#command line option che chiede il csv in entrata

parser = argparse.ArgumentParser(description='Specify the input csvfile and the output excel file', usage='python3 to_csv.py --input file.xlsx --output file.csv')
parser.add_argument('-i', '--input', action='store', type=str, help='-[required] specify the input excel')
parser.add_argument('-o', '--output', action='store', type=str, help='-[required] specify the output csv file name')
args = parser.parse_args()

excelfile=args.input
csvfile=args.output

df = pd.read_excel(excelfile)
df.to_csv(csvfile, sep = ",", index = False)