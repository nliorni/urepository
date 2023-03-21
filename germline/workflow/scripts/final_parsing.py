import pandas as pd
import numpy as np
import argparse

#command line option che chiede il csv in entrata

parser = argparse.ArgumentParser(description='Specify the input csv/tsv file and the output excel file', usage='python3 final_parsing.py --input file.csv --output file.xlsx')
parser.add_argument('-i', '--input', action='store', type=str, help='-[required] specify the input csvfile')
parser.add_argument('-g', '--genelists', nargs = "*", help='-specify input gene lists in .txt format', required = False)
parser.add_argument('-o', '--output', action='store', type=str, help='-[required] specify the output excel file name')
args = parser.parse_args()

csv_input=args.input
gene_list_paths = args.genelists
excel_output=args.output

def get_genotype(df):
    
    if df["HOM"]:
        return "hom"
    elif df["HET"]:
        return "het"
    else:
        return ""

def predictors_check(df, threshold_dict):
    
    # initialize new column
    df["PredictorsCount"] = "."
    
    for index, row in df.iterrows():
        
        # initialize counter
        
        all_predictors = str(len(threshold_dict))
        
        exceeds_threshold = 0
        
        # cycle over the predictors
        for predictor, threshold in threshold_dict.items():
            if isinstance(threshold, list):
                if row[predictor] in threshold:
                    exceeds_threshold += 1
                    #break
            elif isinstance(threshold, str):
                if row[predictor] == threshold:
                    exceeds_threshold += 1
                    #break
            else:
                if row[predictor] > threshold:
                    exceeds_threshold += 1
                    #break
        
        exceeds_threshold = str(exceeds_threshold)
                    
        # store in the new column
        df.at[index, "PredictorsCount"] = f"{exceeds_threshold}/{all_predictors}"
    
    return(df)

# Read Excel file
print("loading file...")
df = pd.read_csv(csv_input, sep = "\t")

# Retrieve headers of the df
header = list(df.columns)

print("checking columns...")
if 'GERP++_NR' in header:
    df = df = df.rename(columns={'GERP++_NR': 'GERP++_NR_'})
if 'GERP++_RS' in header:
    df = df.rename(columns={'GERP++_RS': 'GERP++_RS_'})
if 'Aloft_prob_Dominant' in header:
    df = df.rename(columns={'Aloft_prob_Dominant': 'Aloft_prob_Dominant_'})
if 'Aloft_prob_Recessive' in header:
    df = df.rename(columns={'Aloft_prob_Recessive': 'Aloft_prob_Recessive_'})
if 'Aloft_prob_Tolerant' in header:
    df = df.rename(columns={'Aloft_prob_Tolerant': 'Aloft_prob_Tolerant_'})


info = ['CHROM', 'POS', 'REF', 'ALT', 'QUAL', 'HOM', 'HET', 'DP', 'GEN[*].DP', 'GEN[*].AD', 'GEN[*].GQ', 'GEN[*].PL', 'FS', 'FILTER', 'Existing_variation', 'Consequence', 'SYMBOL', 'Feature_type', 'Feature', 'BIOTYPE', 'EXON', 'INTRON', 'HGVSc', 'HGVSp', 'cDNA_position', 'CDS_position', 'Protein_position', 'Amino_acids', 'Codons', 'STRAND', 'CADD_PHRED', 'ClinVar', 'ClinVar_CLNSIG', 'ClinVar_CLNDN', 'gnomAD_211_AF', 'gnomAD_211_AF_afr', 'SIFT', 'PolyPhen', 'miRNA', 'APPRIS']
aloft = ['Aloft_Confidence', 'Aloft_Fraction_transcripts_affected', 'Aloft_pred', 'Aloft_prob_Dominant_', 'Aloft_prob_Recessive_', 'Aloft_prob_Tolerant_']
multival = ['BayesDel_addAF_pred', 'BayesDel_addAF_score', 'BayesDel_noAF_pred', 'ClinPred_pred', 'ClinPred_score', 'DANN_score', 'DEOGEN2_pred', 'DEOGEN2_score', 'Eigen-PC-phred_coding', 'Eigen-phred_coding', 'FATHMM_pred', 'FATHMM_score', 'GERP++_NR_', 'GERP++_RS_', 'GenoCanyon_score', 'LIST-S2_pred', 'LIST-S2_score', 'LRT_pred', 'LRT_score', 'M-CAP_pred', 'M-CAP_score', 'MPC_score', 'MVP_score', 'MetaLR_pred', 'MetaLR_score', 'MetaSVM_pred', 'MetaSVM_score', 'MutPred_score', 'MutationAssessor_pred', 'MutationAssessor_score', 'MutationTaster_pred', 'MutationTaster_score', 'PROVEAN_pred', 'PROVEAN_score', 'Polyphen2_HDIV_pred', 'Polyphen2_HDIV_score', 'Polyphen2_HVAR_pred', 'Polyphen2_HVAR_score', 'REVEL_score','Reliability_index', 'SIFT4G_pred', 'SIFT4G_score', 'SIFT_pred', 'SIFT_score', 'VEST4_score']
fixed = ['MutPred_Top5features', 'Interpro_domain','SiPhy_29way_logOdds', 'SiPhy_29way_pi', 'fathmm-MKL_coding_pred', 'fathmm-MKL_coding_score', 'fathmm-XF_coding_pred', 'fathmm-XF_coding_score', 'integrated_confidence_value', 'integrated_fitCons_score', 'phastCons100way_vertebrate', 'phastCons17way_primate', 'phastCons30way_mammalian', 'phyloP100way_vertebrate', 'phyloP17way_primate', 'phyloP30way_mammalian', 'ExAC_AF', 'ExAC_AFR_AF', 'ESP6500_AA_AF', 'ESP6500_EA_AF']


print("removing &s from INFOs...")
df_info = df.loc[:,info]
## Selecting first Transcipt (Existing_variation) and the correpsonding APRIS
print("selecting first transcript...")
df_info['Existing_variation'] = df_info['Existing_variation'].str.split('&').str[0].str.split(',').str[0]
df_info['APPRIS'] = df_info['APPRIS'].str.split('&').str[0]
df_info['Consequence'] = df_info['Consequence'].str.split('&').str[0]
df_info['ClinVar_CLNDN'] = df_info['ClinVar_CLNDN'].str.split('&').str[0]

print("removing &s from Aloft...")
df_aloft = df.loc[:,aloft]
header_aloft = list(df_aloft.columns)
df_aloft = df_aloft.replace(",", "&", regex=True)

print("removing &s from multival...")
df_multival = df.loc[:,multival]
header_multival = list(df_multival.columns)
df_multival = df_multival.replace(",", "&", regex=True)
df_multival = df_multival.replace("'", "", regex=True)

print("removing &s from fixed...")
df_fixed = df.loc[:,fixed]
df_fixed['Interpro_domain'] = df_fixed['Interpro_domain'].str.split('&').str[0]

effect = np.array(["High", "A", "D", "Dominant", "Recessive" ,"H", "M", "L", "P", "N", "T", "Tolerant", "B", "U", "Low", ".", "-", ""])

print("some format operations...")
## Seprate in df_multival the coupled columns from the independent ones
coupled_columns = df_multival.columns[df_multival.columns.str.replace("_pred", "").str.replace("_score", "").duplicated(keep=False)].values.reshape(-1,2)
single_columns = df_multival.columns.difference(coupled_columns.reshape(-1))

## Separate the lietral and numerical ones
alfa_col = []
for col in df_multival.columns:
  break_outer_loop = False
  for element in df_multival[col].astype(str).str.split("&").values:
    for x in element:
      if x in effect[:-3]:
        alfa_col.append(col)
        break_outer_loop = True
        break
    if break_outer_loop:
      break

alfa_col = pd.Index(alfa_col)
numeric_col = df_multival.columns.difference(alfa_col)

print("handling max scores and max severity...")

def return_max_scores(series):
  series = series.copy().astype(str).str.split("&")
  return pd.to_numeric(series.dropna().apply(max).reindex(series.index), errors="coerce")

df_multival[single_columns.intersection(numeric_col)] = df_multival[single_columns.intersection(numeric_col)].apply(return_max_scores)

def return_max_severity(series):
  series = series.copy()
  return series.dropna().astype(str).str.split("&").apply(lambda lst: sorted(lst, key=effect.tolist().index)[0]).reindex(series.index)

df_multival[single_columns.intersection(alfa_col)] = df_multival[single_columns.intersection(alfa_col)].apply(return_max_severity)

def return_max_severity_and_score(series_pred, series_score):

  starting_index = series_pred.index
  series_pred, series_score = series_pred.copy().dropna(), series_score.copy().dropna()

  series_pred = series_pred.astype(str).str.split("&")

  max_indices = series_pred.apply(lambda lst: lst.index(sorted(lst, key=effect.tolist().index)[0]))
  series_pred = series_pred.apply(lambda lst: sorted(lst, key=effect.tolist().index)[0])  
  series_score = pd.concat([series_score.astype(str).str.split("&"), max_indices], axis=1, ignore_index=True)\
  .apply(lambda s: s.iloc[0][s.iloc[1]], axis=1)
  
  return series_pred.reindex(starting_index), series_score.reindex(starting_index)

for col1, col2 in coupled_columns:
  try:
    df_multival[col1] ,df_multival[col2] = return_max_severity_and_score(df_multival[col1], df_multival[col2])
  except: #the exepcion occours if the columns have only NaN values
    pass

# ####
#Working on the "Aloft" columns
####

for index,row in df_aloft.iterrows():
    for i in range(-1, len(header_aloft)-1):
        c = 0
        col_curr = ''.join(df_aloft.columns[i].strip().split("_")[:-1])
        col_next = ''.join(df_aloft.columns[i+1].strip().split("_")[:-1])
        col_prev = ''.join(df_aloft.columns[i-1].strip().split("_")[:-1])
        #print(col_curr, col_prev, col_next)
        if row[i] == "NaN":
            row[i] == row[i]
            #print(row[i])
        elif str(row[i]).find("&") < 0:
            row[i] == row[i]
            #print(row[i])
        elif str(row[i]).find("&") >= 1:
            val = list(row[i].strip().split("&"))
            for x in val:
                if x.isalpha()==False:
                    c += 1
                if c < 1 :
                    pred = val
                    sorted_pred = sorted(pred,key=effect.tolist().index)
                    row[i] = sorted_pred[0]
                else:
                    row[i] = pd.Series(row[i].strip().split("&")).replace\
                    ({'.': np.nan, '': np.nan, '-': np.nan})\
                    .sort_values(ascending=False).fillna(".").iloc[0]

if 'GERP++_NR_' in list(df_multival.columns):
    df_multival = df_multival.rename(columns={'GERP++_NR_': 'GERP++_NR'})
if 'GERP++_RS_' in list(df_multival.columns):
    df_multival = df_multival.rename(columns={'GERP++_RS_': 'GERP++_RS'})

if 'Aloft_prob_Dominant_' in list(df_aloft.columns):
    df_aloft = df_aloft.rename(columns={'Aloft_prob_Dominant_': 'Aloft_prob_Dominant'})
if 'Aloft_prob_Recessive_' in list(df_aloft.columns):
    df_aloft = df_aloft.rename(columns={'Aloft_prob_Recessive_': 'Aloft_prob_Recessive'})
if 'Aloft_prob_Tolerant_' in list(df_aloft.columns):
    df_aloft = df_aloft.rename(columns={'Aloft_prob_Tolerant_': 'Aloft_prob_Tolerant'})


# Merging the 4 sub-DataFrame 
print("merging info, aloft, multival and fixed columns...")
df_final = pd.concat([df_info, df_aloft, df_multival, df_fixed], axis=1)
df_final = df_final.replace(",-", "", regex=True)
df_final = df_final.replace("&-", "", regex=True)
df_final = df_final.replace(to_replace='^\.$', value=np.nan, regex=True)

rename_dict = {"GEN[*].GQ":"GenotypeQuality",
               "DP":"TotalDepth", 
               "GEN[*].DP":"GenotypeDepth", 
               "GEN[*].AD":"AlleleDepth",
               "GEN[*].PL":"GenotypesLikelihood",
               "HGVSc":"GeneDetail",
               "BIOTYPE":"GeneFunction",
               "Feature":"RefSeqTranscriptId",
               "SYMBOL":"GeneSymbol",
               "CHROM":"chr",
               "POS":"pos",
               "REF":"ref",
               "ALT":"alt",
               "QUAL":"qual",
               "FILTER":"filter",
               "Feature_type":"FeatureType",
               "EXON":"Exon",
               "INTRON":"Intron",
               "cDNA_position":"cDNAposition",
               "CDS_position":"CDSposition",
               "Protein_position":"ProteinPosition",
               "Amino_acids":"AminoAcids",
               "STRAND": "Strand",
               "Interpro_domain":"InterproDomain",
               "Existing_variation":"dbSNP",
               "ClinVar_CLNSIG":"ClinVarCLNSIG",
               "ClinVar_CLNDN":"ClinVarCLNDN",
              }

## TODO: add DANN score... find the correct threshold
predictors = {"CADD_PHRED":15,  
              "SIFT": "deleterious", 
              "ClinPred_pred":"D",
              "DEOGEN2_pred":"D",
              "FATHMM_pred":"D",
              "LIST-S2_pred":"D",
              "LRT_pred":"D",
              "M-CAP_pred":"D",
              "MetaLR_pred":"D",
              "MetaSVM_pred":"D",
              "MutationAssessor_pred":["H","M"],
              "MutationTaster_pred":["A", "D"],
              "PROVEAN_pred":"D",
              "Polyphen2_HDIV_pred":["D", "P"],
              "Polyphen2_HVAR_pred":["D", "P"],
              "fathmm-MKL_coding_pred":"D",
              "fathmm-XF_coding_pred":"D",
             }

print("renaming columns as needed...")
df_final = df_final.rename(columns = rename_dict)
df_final["CADD_PHRED"] = pd.to_numeric(df_final["CADD_PHRED"], errors = "coerce")

print("checking predictors...")
predictors_check(df_final, predictors)
print("formatting genotype column...")
df_final.insert(7, "Genotype", df_final.apply(get_genotype, axis = 1))
df_final = df_final.drop(columns = ["HOM", "HET"])

print("reordering and keeping only useful columns...")
reorder_keep = ["chr", "pos", "ref", "alt", "qual", "filter", "Genotype", "GenotypeQuality", "TotalDepth", "GenotypeDepth", "AlleleDepth", "GenotypesLikelihood", "FS", "GeneSymbol", "GeneFunction", "Consequence", "RefSeqTranscriptId", "Strand", "cDNAposition", "CDSposition", "ProteinPosition", "Codons", "Exon", "Intron", "GeneDetail", "HGVSp", "AminoAcids", "InterproDomain", "dbSNP", "ClinVarCLNSIG", "ClinVarCLNDN", "GERP++_RS", "phyloP100way_vertebrate", "CADD_PHRED", "SIFT", "DANN_score", "DEOGEN2_pred", "Eigen-phred_coding", "FATHMM_pred", "LRT_pred", "M-CAP_pred", "MetaLR_pred", "MetaSVM_pred",  "MutationAssessor_pred", "MutationTaster_pred", "PROVEAN_pred", "Polyphen2_HDIV_pred", "Polyphen2_HVAR_pred", "REVEL_score", "MutPred_Top5features", "fathmm-MKL_coding_pred", "fathmm-XF_coding_pred", "gnomAD_211_AF", "gnomAD_211_AF_afr", "ExAC_AF", "ExAC_AFR_AF", "ESP6500_AA_AF", "ESP6500_EA_AF", "PredictorsCount"]
df_final = df_final[reorder_keep]

print("checking gene lists...")
if gene_list_paths is None:
    print("WARNING: no gene lists were provided, keeping the file untouched...")
    df_final.to_excel(excel_output, index = None)
else:
    for gene_list_path in gene_list_paths:
        gene_list_name = gene_list_path.rstrip(".txt").rsplit("/")[-1]
        print("-"+gene_list_name)
        gene_list = pd.read_csv(gene_list_path, sep = "\t", names = ["GeneSymbol"])
        gene_list[gene_list_name] = "X"
        df_final = pd.merge(df_final, gene_list, how = "left", on = "GeneSymbol")

## edit to read a file list of columns to keep
tonum = ["GenotypeQuality", "TotalDepth", "GenotypeDepth", "GERP++_RS", "phyloP100way_vertebrate", "DANN_score", "Eigen-phred_coding", "REVEL_score", "gnomAD_211_AF", "gnomAD_211_AF_afr", "ExAC_AF", "ExAC_AFR_AF", "ESP6500_AA_AF", "ESP6500_EA_AF"]

print("converting string columns to numeric...")
for num in tonum:
    print("-"+num)
    df_final[num] = pd.to_numeric(df_final[num], errors = "coerce")



# clean orphan chars (' - -,)
print("cleaning output file from orphan chars: getting rid of ', -,. ,  - and -,  ...") 
#cleaner = {}
cleanexp = {r'^-,.$':'', r'^,.$': '', r'^,\.$': '', r'^-,\.$': '', r'^,$': '',  r'^-,$': '', r'^-$': ''}
df_final = df_final.replace(cleanexp, regex = True)


#Exporting into an output Excel File
print("saving output file...")
df_final.to_excel(excel_output, index=False)