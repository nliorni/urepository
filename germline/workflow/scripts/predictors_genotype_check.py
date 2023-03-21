import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='Specify the input excel and the output excel file', usage='python3 select_value.py --excel file.xlsx --output out_file.xlsx')
parser.add_argument('-e', '--excel', action='store', type=str, help='-[required] specify the input excel file name')
parser.add_argument('-o', '--output', action='store', type=str, help='-[required] specify the output excel file name')
args = parser.parse_args()

path = args.excel
output = args.output

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

## TODO: add DANN score
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

print("reading excel input file")
file = pd.read_excel(path)
print("renaming columns as needed")
file = file.rename(columns = rename_dict)
file["CADD_PHRED"] = pd.to_numeric(file["CADD_PHRED"], errors = "coerce")

print("checking predictors")
predictors_check(file, predictors)
print("formatting genotype column")
file.insert(7, "Genotype", file.apply(get_genotype, axis = 1))
file = file.drop(columns = ["HOM", "HET"])

print("reordering and keeping only useful columns")
reorder_keep = ["chr", "pos", "ref", "alt", "qual", "filter", "Genotype", "GenotypeQuality", "TotalDepth", "GenotypeDepth", "AlleleDepth", "GenotypesLikelihood", "GeneSymbol", "GeneFunction", "Consequence", "RefSeqTranscriptId", "Exon", "Intron", "GeneDetail", "AminoAcids", "InterproDomain", "dbSNP", "ClinVarCLNDN", "GERP++_RS", "phyloP100way_vertebrate", "CADD_PHRED", "SIFT", "DANN_score", "DEOGEN2_pred", "Eigen-phred_coding", "FATHMM_pred", "LRT_pred", "M-CAP_pred", "MetaLR_pred", "MetaSVM_pred",  "MutationAssessor_pred", "MutationTaster_pred", "PROVEAN_pred", "Polyphen2_HDIV_pred", "Polyphen2_HVAR_pred", "REVEL_score", "MutPred_Top5features", "fathmm-MKL_coding_pred", "fathmm-XF_coding_pred", "gnomAD_211_AF", "gnomAD_211_AF_afr", "ExAC_AF", "ExAC_AFR_AF", "ESP6500_AA_AF", "ESP6500_EA_AF", "PredictorsCount"]
file = file[reorder_keep]

file.to_excel(output, index = None)
