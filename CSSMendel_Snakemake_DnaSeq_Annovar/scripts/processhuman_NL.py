#!/usr/bin/env python2.7

"""
TITLE
"""
import os
import sys
import getopt
import re
import pymongo
import itertools
import argparse
import pandas
from pandas.errors import EmptyDataError
# import pandas.core.format
import sys
from subprocess import Popen
from dbsnp_annotator import dbsnp_calculator
from ordered_set import OrderedSet

# pandas.set_option('display.max_rows', 5000)
# pandas.set_option('display.max_columns', 500)
# pandas.set_option('display.width', 1000)

__author__ = "Daniele Capocefalo, Mauro Truglio"
__copyright__ = "Copyright 201X"
__credits__ = [""]
__version__ = "0.0.1"
__maintainer__ = "Mauro Truglio"
__email__ = "bioinformatics@css-mendel.it"
__status__ = "Development"
__date__ = ""
__license__ = u"""
  Copyright (C) 2016-2017  Mauro Truglio <m.truglio@css-mendel.it>
  Viale Regina Margherita 261, 00198 Rome, Italy

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA
  02110-1301 USA
  """

def atoi(text):
    return int(text) if text.isdigit() else text


def natural_keys(text):
    """
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    """
    return [atoi(c) for c in re.split('(\d+)', text)]


def listmover(l, newind, oldind):
    # insert element to ind1
    l.insert(newind, l[oldind])
    del l[oldind+1]
    return l


def processhuman(args, rsidind=-1, indlist=[]): #indlist=[7, 13, 14]):
    """
    this module handles .human files created by parseVCF and modify them in order to add variant Annotation
    """
    # first, clean the file of all the useless column that are no longer needed
    print("Removing useless columns from input...\n")
    # remove useless columns from human file
    humanfile = args.filein
    outdir = args.d
    cleanedhuman = os.path.splitext(humanfile)[0] + "_cleaned.human"


    with open(humanfile, "r") as humanin:
        with open(cleanedhuman, "w") as out:
            for i, line in enumerate(humanin):
                tmp = line.strip().split("\t")
                tmp = [i for j, i in enumerate(tmp) if j not in indlist]
                out.write("\t".join(tmp) + "\n")

    print("Input file was cleaned!\n")
    # remove intermediate files
    #os.remove(humanfile)
    
    outsingle = os.path.join(outdir, os.path.splitext(os.path.basename(cleanedhuman))[0] + "_single.human") # path to the single annotations
    outmulti = os.path.join(outdir, os.path.splitext(os.path.basename(cleanedhuman))[0] + "_multiple_human.tsv")  # path to the multiple annotations
    # print outmulti
    # insert the PATHS of the two output in a string
    outlist = [outsingle, outmulti]
    # print outlist
    coordic = {}
    with open(cleanedhuman,"r") as human:  # open all files for writing since they will be used for these operations accordingly
        with open(outsingle, "w") as singles:
            with open(outmulti, "w") as multi:
                # write header
                header = human.readline()
                singles.write(header)
                multi.write(header)

                for line in human:
                    tmp = line.split("\t")
                    coords = "\t".join(["\t".join(tmp[0:3]), tmp[rsidind]])
                    if coords not in coordic:
                        coordic[coords] = []
                        coordic[coords].append(line)
                    else:
                        coordic[coords].append(line)

                for elem in coordic.keys():
                    if len(coordic[elem]) >= 2:
                        # print coordic[elem]
                        multi.write("".join(coordic[elem]))
                    else:
                        singles.write("".join(coordic[elem]))
    return outlist


def annotate_variants(humanfile, outdir, annoaction, splicing_window, completeflag, annovarcmd, humandb, customfilters=None, buildv="hg19"):
    '''
    This module covers the annotation of variants already parsed using the table.annovar.pl script on a series of references sequences
    '''

    print("Proceeding with ANNOVAR Annotation according to your tastes...\n")
    outfile = os.path.join(outdir, os.path.splitext(os.path.basename(humanfile))[0])
    realout = outfile + ".{}_multianno.txt".format(buildv)  # path of the real output file
    extrafields = "--nastring . --dot2underline"
    print("INFO:", humanfile, annoaction)
    if "cleaned_single" in humanfile:
        tempdir = outdir
    else:
        if os.path.exists("/scratch/tmp"):
            tempdir = "/scratch/tmp"
        else:
            tempdir = "/tmp"
    print("Using", tempdir, "as temp dir")
    stdargs = " ".join([annovarcmd, extrafields, "--out", outfile, "--buildver ", buildv, "--tempdir", tempdir])

    finalshell = ""
    finalfile = ""

    if annoaction == "standard":
        if completeflag == "True":  # print out all fields
            if buildv=="hg19":
                prots = ["refGene", "clinvar_20210123", "spidex", "dbscsnv11", "dbnsfp40aTranscripts", "cadd14", "gnomad211_exome",
                         "esp6500siv2_all", "esp6500siv2_ea", "esp6500siv2_aa", "exac03", "exac03nontcga", "exac03nonpsych", "hrcr1", "gme"]
            elif buildv=="hg38":
                prots = ["refGene", "clinvar_20210123", "dbscsnv11", "dbnsfp40aTranscripts", "cadd15", "gnomad211_exome",
                         "esp6500siv2_all", "esp6500siv2_ea", "esp6500siv2_aa", "exac03", "exac03nontcga", "exac03nonpsych", "hrcr1", "gme"]
        else:  # print out a limited annotation set
            if buildv=="hg19":
                prots = ["refGene", "clinvar_20210123", "spidex", "dbscsnv11", "dbnsfp40aTranscripts", "cadd14", "gnomad211_exome", "esp6500siv2_all", "exac03"]
            elif buildv == "hg38":
                prots = ["refGene", "clinvar_20210123", "dbscsnv11", "dbnsfp40aTranscripts", "cadd15", "gnomad211_exome", "esp6500siv2_all", "exac03"]


        # print prots
        protocols = " ".join(["--protocol", ",".join(prots)])
        # print protocols
        filters = ",".join(["f"] * (len(prots) - 1))
        # print filters
        operations = "--operation g," + filters
        # print operations
        # define final file based on SAMPLE name
        finalfile = os.path.join(outdir, os.path.splitext(os.path.basename(humanfile))[0] + "_annotated_standard.txt")

    elif annoaction == "somatic":
        if completeflag == "True":
            if buildv=="hg19":
                prots = ["refGene", "clinvar_20210123", "spidex", "dbscsnv11", "dbnsfp40aTranscripts", "cadd14", "nci60", "cosmic70", "icgc21", "gnomad211_exome",
                         "esp6500siv2_all", "esp6500siv2_ea", "esp6500siv2_aa", "exac03",
                         "exac03nontcga", "exac03nonpsych", "kaviar_20150923", "hrcr1", "gme"]
            elif buildv == "hg38":
                prots = ["refGene", "clinvar_20210123", "dbscsnv11", "dbnsfp40aTranscripts", "cadd15", "nci60", "cosmic70", "gnomad211_exome",
                         "esp6500siv2_all", "esp6500siv2_ea", "esp6500siv2_aa", "exac03",
                         "exac03nontcga", "exac03nonpsych", "kaviar_20150923", "hrcr1", "gme"]

        else:
            if buildv=="hg19":
                prots = ["refGene", "clinvar_20210123", "spidex", "dbscsnv11", "dbnsfp40aTranscripts", "cadd14", "nci60", "cosmic70", "icgc21", "gnomad211_exome",
                     "esp6500siv2_all", "exac03"]
            elif buildv == "hg38":
                prots = ["refGene", "clinvar_20210123", "dbscsnv11", "dbnsfp40aTranscripts", "cadd15", "nci60", "cosmic70", "gnomad211_exome",
                     "esp6500siv2_all", "exac03"]

        protocols = " ".join(["--protocol", ",".join(prots)])
        # print protocols
        filters = ",".join(["f"] * (len(prots) - 1))
        # print filters
        operations = "--operation g," + filters
        # print operations
    
        # define final file based on SAMPLE name
        finalfile = os.path.join(outdir, os.path.splitext(os.path.basename(humanfile))[0] + "_annotated_somatic.txt")

    elif annoaction == "custom":  # customfilters is a comma separated string of specific filters
        # remove all empty spaces
        prots = customfilters.split(",")
        print(prots)
        # create list from string of protocols
        prots.insert(0, "refGene")  # add refGene to the protocols list
        protocols = " ".join(["--protocol", ",".join(prots)])
        # print protocols
        operations = " ".join(["--operation", "g," + "f," * (len(prots) - 1), " "])
    
        # define final file based on SAMPLE name
        finalfile =os.path.join(outdir, os.path.splitext(os.path.basename(humanfile))[0] + "_annotated_custom.txt")

    else:
        print("ERROR! action of annotation invalid! could not complete annotation!\nQuitting...\n")
        sys.exit(1)

    if splicing_window != "2":
        splstring = "\'--splicing_threshold" + " " + splicing_window + "\'"
        # print (splstring)
        other = list(itertools.repeat("\'\'", (len(prots) - 1)))
        splstring = [splstring] + other
        # print (splstring)
        splstring = ",".join(splstring)
        # print(splstring)
        # print(splstring)
        stdargs = " ".join([stdargs, "--argument", splstring])
        # print(stdargs)
    
        finalshell = " ".join([stdargs, protocols, operations, humanfile, humandb])

    else:
        finalshell = " ".join([stdargs, protocols, operations, humanfile, humandb])

    if finalshell and finalfile:

        print("ANNOVAR COMMAND:")
        print(finalshell)
        os.system(finalshell)

    else:
        print("ERROR! CMD OR FINAL FILE NOT DEFINED!\nQuitting...")
        print("SHELL ERROR: " + finalshell)
        print("FINAL PATH: " + finalfile)

        sys.exit(1)

    try:
        os.rename(realout, finalfile)
        print("Annotation was successful! final annotation file  path is: " + os.path.abspath(finalfile))
    except OSError:
        print("ERROR:! ANNOTATION WAS NOT  POSSIBLE, CHECK SHELL LOG FOR MORE DETAILS!\nending..\n")
        sys.exit(1)

    return finalfile


def annomanager(humanfile, args, multiplicity):
    """
    given a list of annotated variants, divide them and extract the minimum information that can be used for annovar,
    then re-add the information straightaway after  parsing
    """

    #load ENSEMBL to REFSEQ dictionary
    ids_dict = {}
    with open("/software/annovar/humandb/ENS_to_refseq.csv") as f:
        next(f)
        for line in f:
            chr = line.split(',')[1].split('_')[0]
            if chr not in ids_dict:
                ids_dict[chr] = {}

            ids_dict[chr][line.split(',')[0]] = line.split(',')[2].strip()

    vardic = {}  # this will serve for reconstruct the file after annovar query
    outdir = args.d
    annoaction = args.action
    samplename = os.path.splitext(os.path.basename(humanfile))[0]
    print('Samplename: {}'.format(samplename))
    splthres = args.splicing_window
    completeflag = args.completeflag
    annovarcmd = args.annovarcmd
    humandb = args.humandb
    customlist = args.customlist
    fields_order = args.fields_order
    buildv = args.genome_version
    sys.stdout.write("Columns reorder requested: {}.".format(fields_order))
    # single variants only for the moment
    # open input file and divide the header from the rest


    with open(humanfile, "r") as inv:
    
        invar = inv.readlines()

    
        prepath = os.path.join(outdir, os.path.splitext(os.path.basename(humanfile))[0] + "_preanno.human")
        prefile = open(prepath, "w")
        
        # parse variants one by one and save them to a dictionary of variants
        for line in invar[1:]:
            tmp = line.rstrip().split("\t")
            var = "\t".join(tmp[:5])

            # print "\t".join(tmp[5:])

            vardic[var] = "\t".join(tmp[5:])
            prefile.write(var + "\n")

        prefile.close()

        
    # now that the output file contains only the needed variants, let's call annovar
    if annoaction == "standard" or annoaction == "somatic":
        annout = annotate_variants(prepath, outdir, annoaction, splthres, completeflag, annovarcmd, humandb, None, buildv)
        annovar_temp = os.path.join(outdir, next(os.walk(outdir))[1][0], 'temp.refGene.exonic_variant_function')
        annovar_temp_region = os.path.join(outdir, next(os.walk(outdir))[1][0], 'temp.refGene.variant_function')


        # This takes care of the synonymous/nonsynonymous annotation of AA changes
        dup_dict = {}
        try:
            df_annovartemp = pandas.read_csv(annovar_temp, sep='\t', low_memory=False, header=None)
            lines = df_annovartemp[0]
            df_duplicates = df_annovartemp[lines.isin(lines[lines.duplicated()])]


            for i, row in df_duplicates.iterrows():
                line = str(int(df_duplicates.at[i, 0].split('line')[1]) - 1)
                if df_duplicates.at[i, 1] != 'unknown':
                    if line not in dup_dict:
                        dup_dict[line] = {}
                    dup_dict[line][df_duplicates.at[i, 1]] = df_duplicates.at[i, 2].rstrip(',')

        except EmptyDataError:
            print("\nNote: Exonic variant functions file is empty, no annotation will be provided for Exonic Function and AA Change")


        # This takes care of the region/gene annotation
        regions_dict = {}
        with open(annovar_temp_region) as rg:
            for line in rg:
                unique_id = '-'.join(line.split('\t')[2:]).strip()
                if unique_id not in regions_dict:
                    regions_dict[unique_id] = {}
                regions_dict[unique_id][line.split('\t')[0]] = line.split('\t')[1]


        df_big = pandas.read_csv(annout, sep='\t', encoding='utf-8', index_col=None)

        # Creating a column that contains only the region function with highest priority (for easy filtering)

        simplified_region_tempdf = df_big['Func_refGene'].str.split(';', expand=True)
        if not simplified_region_tempdf.dropna().empty:
            df_big.insert(4, 'Main region', simplified_region_tempdf[0])

        for i, row in df_big.iterrows():

            # Assigning a function to each AA change, in order to disambiguate annovar results.
            if ';' in df_big.at[i, 'ExonicFunc_refGene']:
                if str(i) in dup_dict:
                    print("Replacing in position", i, df_big.at[i, 'ExonicFunc_refGene'])

                    df_big.at[i, 'ExonicFunc_refGene'] = df_big.at[i, 'ExonicFunc_refGene'].split(';')[0].replace(' ',
                                                                                                                  '_')
                    new_effects = []
                    for key in dup_dict[str(i)]:
                        new_effects.append("{0}[{1}]".format(key.replace(' ', '_'), dup_dict[str(i)][key]))
                    df_big.at[i, 'AAChange_refGene'] = '; '.join(new_effects)
                else:
                    sys.exit("ERROR! line %s had ; but does not appear in the dup_dict dictionary!\n" % i)
            elif df_big.at[i, 'ExonicFunc_refGene'] != '.' and df_big.at[i, 'ExonicFunc_refGene'] != 'unknown':
                df_big.at[i, 'ExonicFunc_refGene'] = df_big.at[i, 'ExonicFunc_refGene'].replace(' ', '_')
                new_effects = "{0}[{1}]".format(df_big.at[i, 'ExonicFunc_refGene'].replace(' ', '_'),
                                                df_big.at[i, 'AAChange_refGene'])
                df_big.at[i, 'AAChange_refGene'] = new_effects


            #Assigning a region function to each gene name in order to disambiguate annovar results.
            unique_id = df_big.at[i, 'Chr'] + '-' + str(df_big.at[i, 'Start']) + '-' + str(df_big.at[i, 'End']) + '-' + df_big.at[i, 'Ref'] + '-' + df_big.at[i, 'Alt']
            if unique_id in regions_dict:
                # print(unique_id,'is in regions dict!')
                new_region = ''
                for key in regions_dict[unique_id]:
                    new_region += '{0}[{1}];'.format(key, regions_dict[unique_id][key])
                df_big.at[i, 'Func_refGene'] = new_region.rstrip(';')



            ens_ids = df_big.at[i, 'Ensembl_transcriptid'].split(';')
            current_chr = df_big.at[i, 'Chr']
            new_ids = []
            for ens_id in ens_ids:
                if ens_id != '.':
                    try:
                        refseq_id = ids_dict[current_chr][ens_id]
                    except KeyError:
                        pass
                    else:
                        new_ids.append("{0}({1})".format(ens_id, refseq_id))

            df_big.at[i, 'Ensembl_transcriptid'] = '; '.join(new_ids)

        df_humanfreq = pandas.read_csv(humanfile, sep='\t', encoding='utf-8', index_col=None)

        # Renaming and dropping columns from main df
        df_big.drop(['Chr', 'Start', 'End', 'Ref', 'Alt'], inplace=True, axis=1)
        df_big.rename(columns={'Func_refGene': 'Region', 'Gene_refGene': 'Gene Name', 'GeneDetail_refGene': 'Gene Detail', 'ExonicFunc_refGene': 'Exonic Function',
                               'AAChange_refGene':'AA Change'}, inplace=True)
        
        df_humanfreq.rename(columns={'Total depth': 'Total Depth (Filtered)'}, inplace=True)
        df_humanfreq.rename(columns={'Genotype depth': 'Genotype Depth (Filtered)'}, inplace=True)

        # joining the previously created human+frequencies file with the main df
        df_joined = pandas.concat([df_humanfreq, df_big], axis=1)


        # Default columns reordering step
        cols = list(df_joined)
        AAchange_idx = cols.index('AA Change')
        cols.insert(AAchange_idx+1, cols.pop(cols.index('Interpro_domain')))
        CLINSIG_idx = cols.index('CLNSIG')
        cols.insert(CLINSIG_idx-1, cols.pop(cols.index('dbSNP ID')))
        cols.insert(CLINSIG_idx-1, cols.pop(cols.index('Freq Ref (dbSNP{0})'.format(args.dbsnp_ver))))
        cols.insert(CLINSIG_idx-1, cols.pop(cols.index('Freq Alt (dbSNP{0})'.format(args.dbsnp_ver))))
        df_joined = df_joined[cols]
        # sys.stdout.write("######### COMPLETEFLAG: "+completeflag)
        # print(df_joined.columns)
        # print(len(df_joined.columns))

        # Custom reordering (if --fields-order list was provided):
        if fields_order != 'None':
            fields_order = str(args.fields_order.strip("'").strip('"'))
            old_order = list(df_joined)
            # print("ORIGINAL")
            # print(old_order)
            # print("REQUEST IS")
            # print(fields_order)
            # print("splitted")
            # print(fields_order.split(";"))

            print(list(OrderedSet(fields_order.split(';'))))
            print("RESULT")
            new_order = list(OrderedSet(fields_order.split(';')) & OrderedSet(list(df_joined)))
            print(new_order)
            if new_order == []:
                new_order = old_order
                sys.stdout.write("The new order of columns that you specified does not contain any of the actual columns"
                                 "in the output files. Will restore the default order.")
            df_joined = df_joined[new_order]
            sys.stdout.write("Columns reordered.")
        else:
            sys.stdout.write("No columns reordering performed.")

        if completeflag == "False":
            for c in df_joined.columns:
                if 'rankscore' in c:
                    df_joined.drop(c, inplace=True, axis=1)
            for c in ['bStatistic_rankscore', 'codon_degeneracy', 'codonpos', 'DANN_rankscore', 'DEOGEN2_rankscore', 'Eigen-PC-raw_coding_rankscore', 'Eigen-raw_coding_rankscore', 'ESP6500_AA_AC', 'ESP6500_EA_AC', 'ExAC_AC', 'ExAC_Adj_AC', 'ExAC_Adj_AF', 'ExAC_AFR_AC', 'ExAC_AFR_AF', 'ExAC_AMR_AC', 'ExAC_AMR_AF', 'ExAC_EAS_AC', 'ExAC_EAS_AF', 'ExAC_FIN_AC', 'ExAC_FIN_AF', 'ExAC_NFE_AC', 'ExAC_NFE_AF', 'ExAC_nonpsych_AC', 'ExAC_nonpsych_Adj_AC', 'ExAC_nonpsych_Adj_AF', 'ExAC_nonpsych_AF', 'ExAC_nonpsych_AFR_AC', 'ExAC_nonpsych_AFR_AF', 'ExAC_nonpsych_AMR_AC', 'ExAC_nonpsych_AMR_AF', 'ExAC_nonpsych_EAS_AC', 'ExAC_nonpsych_EAS_AF', 'ExAC_nonpsych_FIN_AC', 'ExAC_nonpsych_FIN_AF', 'ExAC_nonpsych_NFE_AC', 'ExAC_nonpsych_NFE_AF', 'ExAC_nonpsych_SAS_AC', 'ExAC_nonpsych_SAS_AF', 'ExAC_nonTCGA_AC', 'ExAC_nonTCGA_Adj_AC', 'ExAC_nonTCGA_Adj_AF', 'ExAC_nonTCGA_AFR_AC', 'ExAC_nonTCGA_AFR_AF', 'ExAC_nonTCGA_AMR_AC', 'ExAC_nonTCGA_AMR_AF', 'ExAC_nonTCGA_EAS_AC', 'ExAC_nonTCGA_EAS_AF', 'ExAC_nonTCGA_FIN_AC', 'ExAC_nonTCGA_FIN_AF', 'ExAC_nonTCGA_NFE_AC', 'ExAC_nonTCGA_NFE_AF', 'ExAC_nonTCGA_SAS_AC', 'ExAC_nonTCGA_SAS_AF', 'ExAC_SAS_AC', 'ExAC_SAS_AF', 'FATHMM_converted_rankscore', 'fathmm-MKL_coding_rankscore', 'fathmm-XF_coding_rankscore', 'GenoCanyon_rankscore', 'GERP++_NR', 'GERP++_RS_rankscore', 'gnomAD_exomes_AFR_AC', 'gnomAD_exomes_AFR_AF', 'gnomAD_exomes_AFR_AN', 'gnomAD_exomes_AFR_nhomalt', 'gnomAD_exomes_AMR_AC', 'gnomAD_exomes_AMR_AF', 'gnomAD_exomes_AMR_AN', 'gnomAD_exomes_AMR_nhomalt', 'gnomAD_exomes_AN', 'gnomAD_exomes_ASJ_AC', 'gnomAD_exomes_ASJ_AF', 'gnomAD_exomes_ASJ_AN', 'gnomAD_exomes_ASJ_nhomalt', 'gnomAD_exomes_controls_AFR_AC', 'gnomAD_exomes_controls_AFR_AF', 'gnomAD_exomes_controls_AFR_AN', 'gnomAD_exomes_controls_AFR_nhomalt', 'gnomAD_exomes_controls_AMR_AC', 'gnomAD_exomes_controls_AMR_AF', 'gnomAD_exomes_controls_AMR_AN', 'gnomAD_exomes_controls_AMR_nhomalt', 'gnomAD_exomes_controls_AN', 'gnomAD_exomes_controls_ASJ_AC', 'gnomAD_exomes_controls_ASJ_AF', 'gnomAD_exomes_controls_ASJ_AN', 'gnomAD_exomes_controls_ASJ_nhomalt', 'gnomAD_exomes_controls_EAS_AC', 'gnomAD_exomes_controls_EAS_AF', 'gnomAD_exomes_controls_EAS_AN', 'gnomAD_exomes_controls_EAS_nhomalt', 'gnomAD_exomes_controls_FIN_AC', 'gnomAD_exomes_controls_FIN_AF', 'gnomAD_exomes_controls_FIN_AN', 'gnomAD_exomes_controls_FIN_nhomalt', 'gnomAD_exomes_controls_NFE_AC', 'gnomAD_exomes_controls_NFE_AF', 'gnomAD_exomes_controls_NFE_AN', 'gnomAD_exomes_controls_NFE_nhomalt', 'gnomAD_exomes_controls_nhomalt', 'gnomAD_exomes_controls_POPMAX_AC', 'gnomAD_exomes_controls_POPMAX_AF', 'gnomAD_exomes_controls_POPMAX_AN', 'gnomAD_exomes_controls_POPMAX_nhomalt', 'gnomAD_exomes_controls_SAS_AC', 'gnomAD_exomes_controls_SAS_AF', 'gnomAD_exomes_controls_SAS_AN', 'gnomAD_exomes_controls_SAS_nhomalt', 'gnomAD_exomes_EAS_AC', 'gnomAD_exomes_EAS_AF', 'gnomAD_exomes_EAS_AN', 'gnomAD_exomes_EAS_nhomalt', 'gnomAD_exomes_FIN_AC', 'gnomAD_exomes_FIN_AF', 'gnomAD_exomes_FIN_AN', 'gnomAD_exomes_FIN_nhomalt', 'gnomAD_exomes_flag', 'gnomAD_exomes_NFE_AC', 'gnomAD_exomes_NFE_AF', 'gnomAD_exomes_NFE_AN', 'gnomAD_exomes_NFE_nhomalt', 'gnomAD_exomes_nhomalt', 'gnomAD_exomes_POPMAX_AC', 'gnomAD_exomes_POPMAX_AF', 'gnomAD_exomes_POPMAX_AN', 'gnomAD_exomes_POPMAX_nhomalt', 'gnomAD_exomes_SAS_AC', 'gnomAD_exomes_SAS_AF', 'gnomAD_exomes_SAS_AN', 'gnomAD_exomes_SAS_nhomalt', 'gnomAD_genomes_AFR_AC', 'gnomAD_genomes_AFR_AF', 'gnomAD_genomes_AFR_AN', 'gnomAD_genomes_AFR_nhomalt', 'gnomAD_genomes_AMR_AC', 'gnomAD_genomes_AMR_AF', 'gnomAD_genomes_AMR_AN', 'gnomAD_genomes_AMR_nhomalt', 'gnomAD_genomes_AN', 'gnomAD_genomes_ASJ_AC', 'gnomAD_genomes_ASJ_AF', 'gnomAD_genomes_ASJ_AN', 'gnomAD_genomes_ASJ_nhomalt', 'gnomAD_genomes_controls_AC', 'gnomAD_genomes_controls_AF', 'gnomAD_genomes_controls_AFR_AC', 'gnomAD_genomes_controls_AFR_AF', 'gnomAD_genomes_controls_AFR_AN', 'gnomAD_genomes_controls_AFR_nhomalt', 'gnomAD_genomes_controls_AMR_AC', 'gnomAD_genomes_controls_AMR_AF', 'gnomAD_genomes_controls_AMR_AN', 'gnomAD_genomes_controls_AMR_nhomalt', 'gnomAD_genomes_controls_AN', 'gnomAD_genomes_controls_ASJ_AC', 'gnomAD_genomes_controls_ASJ_AF', 'gnomAD_genomes_controls_ASJ_AN', 'gnomAD_genomes_controls_ASJ_nhomalt', 'gnomAD_genomes_controls_EAS_AC', 'gnomAD_genomes_controls_EAS_AF', 'gnomAD_genomes_controls_EAS_AN', 'gnomAD_genomes_controls_EAS_nhomalt', 'gnomAD_genomes_controls_FIN_AC', 'gnomAD_genomes_controls_FIN_AF', 'gnomAD_genomes_controls_FIN_AN', 'gnomAD_genomes_controls_FIN_nhomalt', 'gnomAD_genomes_controls_NFE_AC', 'gnomAD_genomes_controls_NFE_AF', 'gnomAD_genomes_controls_NFE_AN', 'gnomAD_genomes_controls_NFE_nhomalt', 'gnomAD_genomes_controls_nhomalt', 'gnomAD_genomes_controls_POPMAX_AC', 'gnomAD_genomes_controls_POPMAX_AF', 'gnomAD_genomes_controls_POPMAX_AN', 'gnomAD_genomes_controls_POPMAX_nhomalt', 'gnomAD_genomes_EAS_AC', 'gnomAD_genomes_EAS_AF', 'gnomAD_genomes_EAS_AN', 'gnomAD_genomes_EAS_nhomalt', 'gnomAD_genomes_FIN_AC', 'gnomAD_genomes_FIN_AF', 'gnomAD_genomes_FIN_AN', 'gnomAD_genomes_FIN_nhomalt', 'gnomAD_genomes_flag', 'gnomAD_genomes_NFE_AC', 'gnomAD_genomes_NFE_AF', 'gnomAD_genomes_NFE_AN', 'gnomAD_genomes_NFE_nhomalt', 'gnomAD_genomes_nhomalt', 'gnomAD_genomes_POPMAX_AC', 'gnomAD_genomes_POPMAX_AF', 'gnomAD_genomes_POPMAX_AN', 'gnomAD_genomes_POPMAX_nhomalt', 'integrated_fitCons_rankscore', 'LINSIGHT_rankscore', 'LRT_converted_rankscore', 'LRT_Omega', 'M-CAP_rankscore', 'MetaLR_rankscore', 'MetaSVM_rankscore', 'MPC_rankscore', 'MutationAssessor_rankscore', 'MutationTaster_converted_rankscore', 'MutPred_rankscore', 'MVP_rankscore', 'phastCons100way_vertebrate_rankscore', 'phastCons17way_primate_rankscore', 'phastCons30way_mammalian_rankscore', 'phyloP100way_vertebrate_rankscore', 'phyloP17way_primate_rankscore', 'phyloP30way_mammalian_rankscore', 'Polyphen2_HDIV_rankscore', 'Polyphen2_HVAR_rankscore', 'PrimateAI_rankscore', 'PROVEAN_converted_rankscore', 'SIFT_converted_rankscore', 'SIFT4G_converted_rankscore', 'SIFT4G_pred', 'SIFT4G_score', 'SiPhy_29way_logOdds_rankscore', 'VEST4_rankscore']:
                try:
                    df_joined.drop(c, inplace=True, axis=1)
                except:
                    pass

        out = os.path.join(outdir, samplename + "_{0}_{1}_annotated_final.tsv".format(annoaction, multiplicity))

        df_joined.to_csv(out, sep='\t', encoding='utf-8', index=False)
        
        # pandas.core.format.header_style = None
        writer = pandas.ExcelWriter(out.replace('.tsv', '.xlsx'), engine='xlsxwriter')
        
        #THIS to avoid large file error
        writer.book.use_zip64()
        
        df_joined.to_excel(writer, encoding='utf-8', index=False, sheet_name='report')
        workbook = writer.book
        worksheet = writer.sheets['report']
        worksheet.set_zoom(90)
        worksheet.freeze_panes(1, 0)
        worksheet.autofilter(0, 0, 0, 999)
        header_fmt = workbook.add_format({'align': 'left', 'bold': True, 'bottom': 1, 'left': 1, 'right': 1})
        worksheet.set_row(0, None, header_fmt)
        writer.save()
        
    elif annoaction == "custom":
        sys.stdout.write("CUSTOM SELECTED\n")
        sys.stdout.write("Filters:\n")
        sys.stdout.write(customlist)
        out = None
        # We need to implement what happens with this custom option. The filters are received as a comma separated
        # list correctly, now it's all about using them to retain the right columns only.

    else:
        sys.stdout.write("ERROR! Performed action is not specified, or customlist of annotation was not provided, "
              "please see --help to visualize possible annotation options.\nEnding...\n")
        sys.exit(1)
    
    return out


def print_help():
    print("""
Usage: processhuman.py -f [file] -d [path] [standard | somatic]

 -f : input file
 -d : output directory
    """)
    sys.exit()
    
    
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "filein", help="file you are processing\n")
    parser.add_argument("-d", "--dir", default=os.path.abspath("./"), dest="d",
                        help="directory to be used to dump analysis results\n")
    # add flag for somatic variation (annovar annotatio using NCI and COSMIC),
    # keep all files and to add bam files comig from SOLID
    parser.add_argument("-a", "--action",
                        help="choice the annotation you want from ANNOVAR (standard, somatic, custom). If not specified, annotation will be standard.",
                        choices=['standard', 'somatic', 'custom'], default='standard')
    parser.add_argument("--splicing-window", "-sw", default="2", help="define a splicing window for variant calling if necessary")
    parser.add_argument("-k", "--keep", action="store_true", dest="k", help="Keep all files used for analysis (will be in the /results/ folder\n")
    parser.add_argument("--dbsnp_ver", dest="dbsnp_ver", help="Version of dbsnp to choose in mongo\n")
    parser.add_argument("--diagnosis", action="store_true",
                        help="apply diagnosis settings\n")
    parser.add_argument("--customlist", default=None)
    parser.add_argument("--complete", dest='completeflag', choices=["True", "False"], default="False",
                        help="use this option to print out all the columns for standard and somatic annotype, otherwise several fields will be deleted")
    parser.add_argument("--annovarcmd")
    parser.add_argument("--humandb")
    parser.add_argument("--port", type=int)
    parser.add_argument("--ip")
    parser.add_argument("--coverage", action="store_true",
                        help="find the under-covered segments (under 30) and the percentage of the total gene length covered")
    parser.add_argument("--fields-order", default='None', help="List containing the order of the columns in the output file. Defaults to annovar predefined order.")
    parser.add_argument("--genome-version", default='hg19', help="Human genome version in use")

    return parser.parse_args()


if __name__ == '__main__':
    
    args = parse_args()    #
    # Elaborazione file human con rimozione colonne "inutili"
    humans = processhuman(args)
    print(humans)

    # Humanfreqs sono i file human con aggiunta di frequenze dbSNP
    humanfreqs = dbsnp_calculator(humans[0], args)
    if humanfreqs == -1:
        raise BrokenPipeError("Error while adding dbsnp IDs (monoallelic). Check log.")

    humanfreqs_multi = dbsnp_calculator(humans[1], args)
    if humanfreqs_multi == -1:
        raise BrokenPipeError("Error while adding dbsnp IDs (multi alleles). Check log.")

    print(humanfreqs)
    print("Analysing monoallelic")
    annoout = annomanager(humanfreqs, args, 'monoallelic')
    print("Analysing multiallelic")
    annoout_multi = annomanager(humanfreqs_multi, args, 'multiple')

    #os.popen("rm -r {0}/*.txt {0}/*.human {0}/*_human.tsv" .format(args.d))
    #os.popen("rm -r {0}".format(os.path.join(args.d, next(os.walk(args.d))[1][0])))
    #os.popen("rm -r {0}/*.human {0}/*_human.tsv" .format(args.d))

