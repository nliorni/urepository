import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='Specify the input excel, desired gene lists and the output excel file', usage='python3 gene_list.py -e excel.xlsx -g genelist1.txt ...  --output file.xlsx')
parser.add_argument('-g', '--genelists', nargs = "+", help='-specify input gene lists in .txt format')
parser.add_argument('-e', '--excel', action='store', help='-specify input excel')
parser.add_argument('-o', '--output', action='store', help='-specify output excel')
args = parser.parse_args()


excel_path = args.excel
gene_list_paths = args.genelists
output_path = args.output

print("loading file...")
excel = pd.read_excel(excel_path)

print("checking gene lists...")
if not gene_list_paths:
    print("WARNING: no gene lists were provided")
    excel.to_excel(output_path, index = None)
else:
    for gene_list_path in gene_list_paths:
        gene_list_name = gene_list_path.rstrip(".txt").rsplit("/")[-1]
        print(gene_list_name)
        gene_list = pd.read_csv(gene_list_path, sep = "\t", names = ["GeneSymbol"])
        gene_list[gene_list_name] = "X"
        excel = pd.merge(excel, gene_list, how = "left", on = "GeneSymbol")
    excel.to_excel(output_path, index = None)