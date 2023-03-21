import pandas as pd
import io
import argparse

parser = argparse.ArgumentParser(description='Specify the input excel, cohort vcf and the output excel file', usage='python3 to_excel.py --input file.csv --output file.xlsx')
parser.add_argument('-v', '--vcf', action='store', help='-specify input vcf')
parser.add_argument('-e', '--excel', action='store', help='-specify input excel')
parser.add_argument('-o', '--output', action='store', help='-specify output excel')
args = parser.parse_args()

def read_vcf(path):
    with open(path, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    
    df=pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str},
        sep='\t'
        ).rename(columns={'#CHROM': 'CHROM'})
    return df

cohort_excel_path = args.excel
cohort_vcf_path = args.vcf
output = args.output

print("loading files...")
cohort_excel = pd.read_excel(cohort_excel_path)
cohort_vcf = read_vcf(cohort_vcf_path)

print("storing FORMAT and cohort's members genotype columns...")
colnames = list(cohort_vcf.iloc[:, 8:].keys())
cohort_excel[colnames] = cohort_vcf.iloc[:, 8:]

for col in colnames:
    print(f"storing -{col}")

print("reordering as wished...")
cols_to_move = cohort_excel.iloc[:, 56:]
cohort_excel = cohort_excel.drop(cohort_excel.columns[56:], axis = 1)
final_cohort_excel = pd.concat([cohort_excel.iloc[:, :6], cols_to_move, cohort_excel.iloc[:, 6:]], axis = 1)

final_cohort_excel.to_excel(output, index = None)