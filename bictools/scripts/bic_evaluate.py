# -*- coding: utf-8 -*-
"""
Created on Wed Oct 29 11:59:16 2014

@author: sven
"""

import argparse
import textwrap
import pandas as pd
import glob
import sys
import os
import HTSeq
import re
import numpy as np
import bictools as BICl
import bictools as BICp
import pyopenms as oms


#==============================================================================
# Read and write
#==============================================================================
def read_biceps_log(in_file):
    """
    Reads the biceps mutation log and creates three dictionaries.

    Paramters;
    ----------------------
    in_file: str,
             file that contains the biceps log
             must have the form:
             uniprot_id,tag,orig_sequence,sequence
             Q9H1U9,A162T,ALKCHGIGEYYR,TLKCHGIGEYYR

    Returns:
    ----------------------------
    seq_dic: dic,
            <pep_seq>:<orig_seq> (Example: 'IYIDLCNIFPPDLVR':{'IYINLCNIFPPDLVR': ''})
    bic_mut_dic: dic,
             <pep_seq>:{<uniprot_id>:<mutation_tag>} (Example: IYIDLCNIFPPDLVR: {'Q9C0D7': ['N843D']})
    bic_mut_tag_dic: dic,
            <uniprot_id>:[mutatioN-tag] (Example: 'Q9Y6Y1': ['I978T'])

    """
    df_bic = pd.read_csv(in_file, sep=",")
    bic_dic = {}
    bic_dic_SNP_tags = {}
    filter_id = []
    seq_dic = {}
    k = 0
    for i in df_bic.iterrows():
        #        if i[1]["uniprot_id"] == "Q13885":
        #            break
        if len(i[1]["tag"].split(";")) == 1:
            filter_id.append(0)
            pass
        else:
            #ignore double mutations such as:
            # ['N103D', 'D110N']
            tags = i[1]["tag"].split(";")
            if (tags[0][0] == tags[1][-1]) and (tags[0][-1] == tags[1][0]):
                filter_id.append(1)
                continue
            else:
                filter_id.append(0)

        # dictionary for OMS comparison
        if i[1]["sequence"] in bic_dic:
            if i[1]["uniprot_id"] in bic_dic[i[1]["sequence"]]:
                bic_dic[i[1]["sequence"]][i[1]["uniprot_id"]].append(i[1]["tag"])
            else:
                bic_dic[i[1]["sequence"]][i[1]["uniprot_id"]] = [i[1]["tag"]]
        else:
            bic_dic[i[1]["sequence"]] = {i[1]["uniprot_id"]:[i[1]["tag"]]}

        #sequence dictionary

        if i[1]["sequence"] in seq_dic:
            if i[1]["orig_sequence"] in seq_dic[i[1]["sequence"]]:
                pass
            else:
                seq_dic[i[1]["sequence"]][i[1]["orig_sequence"]] = ""
                k += 1
        else:
            seq_dic[i[1]["sequence"]] = {i[1]["orig_sequence"]:""}

        # dictionary for variationDB comparison
        if i[1]["uniprot_id"] in bic_dic_SNP_tags:
            bic_dic_SNP_tags[i[1]["uniprot_id"]].append(i[1]["tag"])
        else:
            bic_dic_SNP_tags[i[1]["uniprot_id"]] = [i[1]["tag"]]
    df_bic["filter"] = filter_id
    #df_bic[df_bic["filter"] == 0]
    return(seq_dic, bic_dic, bic_dic_SNP_tags)

def CanProVar_to_table(in_file, out_folder):
    """
    Writes the CanProvar results from the fasta DB to a table
    and returns the indexed pandas dataframe
    """
    fasta = HTSeq.FastaReader(in_file)
    ensemble_id = []
    dbsnp_ids = []
    FROM = []
    TO = []
    POS = []
    ID = []
    native_id = []

    #iterate over the fasta file (CanProvar Format)
    #and get the mutations that are written to the description
    #create mutation tags and store them in a dataframe
    for seq in fasta:
        #multiple mutations are seperated by ; in the source file
        split = seq.descr.split(";")
        if split[0] != '':
            # for all mutations generate the specific tag
            # i.e. FROM, TO, POS = A,D, 20
            # keep track of the identifier etc...
            for mutation in split:
                single_mut = mutation.split(":")
                pos = int(re.search("(\d+)", single_mut[1]).groups()[0])
                muts = re.search("([A-Z*-]+)\d+([A-Z*-]+)", single_mut[1]).groups()
                POS.append(pos)
                FROM.append(muts[0])
                TO.append(muts[1])
                ensemble_id.append(seq.name)
                dbsnp_ids.append(single_mut[0])
                native_id.append(single_mut[1])
                ID.append("CanProVar")

    # convert lists to dataframe
    canprovar_df = pd.DataFrame()
    canprovar_df["FROM"] = FROM
    canprovar_df["TO"] = TO
    canprovar_df["POS"] = POS
    canprovar_df["ID"] = ID
    canprovar_df["native_id"] = native_id
    canprovar_df["ensemble_id"] = ensemble_id
    canprovar_df["dbsnp_ids"] = dbsnp_ids

    # extract only the ensemble_ids
    ensemble_ids = pd.DataFrame()
    ensemble_ids["ids"] = np.unique(canprovar_df["ensemble_id"])
    ensemble_ids.to_csv(out_folder+"canprovar_ensemble_ids.csv", sep="\t")

    #index dataframe to be adressed by the ensemble id via
    #canprovar_df.loc["ENSP00000370532"]
    canprovar_df.to_csv(out_folder+"canprovar_tab.csv", sep="\t")
    canprovar_df = canprovar_df.set_index("ensemble_id")
    return(canprovar_df)


def create_mutation_tag(uniprot_id, orig_sequence, sequence, uniprot_dic):
    """
    Creates a mutation tag given a protein identifer and the identified
    peptides. Obviously, this is thought to identify shared peptides
    in different proteins (shared peptides have different positions)
    and thus differnt mutation tags

    Paramters:
    --------------------------------
    """
    #uniprot_id, orig_sequence, sequence, uniprot_dic = item, seq_dic[seqi],
    #seqi, uniprot_dic

    uni_sequence = uniprot_dic[uniprot_id]

    got_sequence = False
    for possible_peptide in orig_sequence.keys():
        if possible_peptide in uni_sequence:
            pos = uni_sequence.index(possible_peptide)
            got_sequence = True
        else:
            continue

        TAG = []
        for i,j,ii in zip(possible_peptide, sequence, range(0, len(sequence) + 1)):
            if i != j:
                TAG.append("%s%s%s" % (i, ii+pos+1, j))

    if got_sequence is False:
        print "error! peptide sequence not found in complete fasta sequence!"
        print "fasta:", orig_sequence.keys()
        print "peptides:", uniprot_dic[item]
        print "protein", item
    return(";".join(TAG))


def adjust_mutation_position(mutation_tags):
    """
    Adjust mutation tags for their identifier as in CanProVar (1-based)
    Add on the biceps based mutation tag +1 as position and we're good

    Paramters:
    --------------------
    mutation-tag: str,
                  tag describing a AA mutation
    """
    #mutation_tags = bic_mut_tag_dic[uniprot_idi]
    new_tags = []
    for item in mutation_tags:
        #double mutations
        items = item.split(";")
        for itemsi in items:
            FROM, POS, TO = re.search("([A-Z*-]+)(\d+)([A-Z*-]+)", itemsi).groups()
            new_tags.append(FROM+str(int(POS)+1)+TO)
    return(new_tags)


def write_to_provean_format(df_val, out_file, mass_filter=1.):
    """
    Given a datafram of validated SAAVs create a valid input list
    for the webservice of PROVEA
    (http://provean.jcvi.org/protein_batch_submit.php?species=human).

    The mass_filter can be used as arbitrary filter for not allowing any
    mutations that might be detected because of low resolution data.

    Paramters:
    ----------------------------------
    df: dataframe,
        contains at least the following columsns ["identifier_uni","tag"]
    """
    fobj = open(out_file, "w")
    fobj.write("uniprot_id POS FROM TO\r\n")
    mass_diffs = []
    for row in df_val.iterrows():
        for tagi in row[1]["tag"].split(";"):
            muts = re.search("([A-Z*-]+)(\d+)([A-Z*-]+)", tagi).groups()
            mdiff = BICp.get_mass_diff(muts[0], muts[2])
            if mdiff >= mass_filter:
                fobj.write("%s %s %s %s\r\n" % (row[1]["identifier_uni"], muts[1], muts[0], muts[2]))
                mass_diffs.append(mdiff)
    fobj.close()


def uniprot_mapping(in_file):
    """
    Reads the uniprot mapping and returns a dictionary with  Uniprot ->ENSEMBLE
    The in_file is retrieved from the uniprot mapping service
    using 4 columns

    Entry\tEntry name\tStatus\tensemble
    Q012T   BLA_HUMAN    reviewed   ENSP0000012

    Returns:
    ----------------------------
    mapping_dic: dic,
            Ensemble uiprot dictionary mapping
    """
    df = pd.read_csv(in_file, sep="\t")
    mapping_dictionary = {}
    for i in df.iterrows():
        #        if i[1]["ensemble"] in mapping_dictionary:
        #            mapping_dictionary[i[1]["ensemble"]].append(i[1]["Entry"])
        #        else:
        #            mapping_dictionary[i[1]["ensemble"]] = [i[1]["Entry"]]
        if i[1]["Entry"] in mapping_dictionary:
            mapping_dictionary[i[1]["Entry"]].append(i[1]["ensemble"])
        else:
            mapping_dictionary[i[1]["Entry"]] = [i[1]["ensemble"]]
    return(mapping_dictionary)


#==============================================================================
# MISC
#==============================================================================
def path_or_dir(in_arg):
    """
    Format arguments depending on it's nature (either path or file)
    Paramter:
    ---------------
    in_arg: str,
            path or file

    Returns:
    ----------------
    in_arg_list: list,
                 list version of the input argument
    """
    if os.path.isfile(in_arg):
        print ("single file mode:")
        in_arg_list = [in_arg]
    elif os.path.isdir(in_arg):
        in_arg_list = sorted(glob.glob(in_arg+"Biceps.Results.*.txt"))
        print ("batch mode: %s files" % len(in_arg_list))
    else:
        print "Error! -in_bic is not a correct file or directory!"
        sys.exit()
    return(in_arg_list)


def main():
    """
    Main function
    """
    parser = argparse.ArgumentParser(
         prog='bic_CSV.py',
         formatter_class=argparse.RawDescriptionHelpFormatter,
         description=textwrap.dedent('''\
         ---------------------------------------------------------------------
         | Description:                                                      |
         ---------------------------------------------------------------------
         Command line tool to convert BICEPS output to a .csv like format.
         The default seperator are tabs ("\\t").

         Explanations for all columns can be found on the github page
         (http://buotex.github.io/BICEPS/doc/html/index.html).

         ---------------------------------------------------------------------
         | Usage:                                                            |
         ---------------------------------------------------------------------
         python BICEPS_TO_CSV.py --infile input.txt --outfile output.csv


         =====================================================================

                 '''))

    parser.add_argument('--in_oms', metavar='in_oms', type=str,
                       help='OpenMS output directory of the TextExporter node.', required=True)
    parser.add_argument('--in_bic', metavar='in_bic', type=str,
                       help='Input file/directory of the processed biceps results  (.mut files)', required=True)
    parser.add_argument('--in_db', metavar='in_db', type=str,
                       help='protein database used for the --in_bic file', required=True)

    parser.add_argument('--CPV_tab', metavar='CPV_tab', type=str,
                       help='Table derived from CanProVar.', default="")
    parser.add_argument('--map_tab', metavar='map_tab', type=str,
                       help='Ensemble-Uniprot mapping data', default="")
    parser.add_argument('--out_aux', metavar='out_aux', type=str,
                       help='Outputfile for auxillary data files', default="")
    parser.add_argument('--min_diff', metavar='min_diff', type=float,
                       help='minimum mass difference between a AA mutation', default=0.5)

    args = parser.parse_args()
    in_oms = args.in_oms
    in_bic = args.in_bic
    in_db = args.in_db
    CPV_tab = args.CPV_tab
    map_tab = args.map_tab
    out_ax = args.out_aux
    min_diff = args.min_diff

    print ""
    print "###################################################################"
    print "# in_oms: %s" % in_oms
    print "# in_bic: %s" % in_bic
    print "# in_db: %s" % in_db
    print "#"
    print "# out dir: %s" % out_ax
    print "###################################################################"
    print ""

    if not os.path.exists(out_ax):
        print "%s created!" % out_ax
        os.makedirs(out_ax)

    #==========================================================================
    #  start analysis workflow
    #==========================================================================
#        in_oms = "/home/sven/data/BICEPS/HeLa/OMS_results2/TOPPAS_out/011-TextExporter-out/"
#        in_bic = "/home/sven/data/BICEPS/HeLa/results/"
#        in_db = '/home/sven/data/BICEPS/fasta/HS_reviewed_19092014.fasta'
#        CPV_tab = "/home/sven/data/BICEPS/var_dbs/Ensembl54_homo_cancer_dbSNP_variation_protein.fasta"
#        map_tab = "/home/sven/data/BICEPS/var_dbs/ensemble_uniprot_mapping_canprovar.tab"
#        out_ax = "/home/sven/data/BICEPS/HCT/results/HCT1/out/"

    print "load fasta files...",
    canprovar_df = CanProVar_to_table(CPV_tab, out_ax)
    mapping_dictionary = uniprot_mapping(map_tab)
    uniprot_dic = BICl.read_uniprot_to_dic(in_db, mode="seq")
    print "DONE!"
    oms_files = sorted(glob.glob(in_oms+"SAAV_candidates*.csv"))
    bic_files = sorted(glob.glob(in_bic+"*.txt.mut"))

    fdr_c = []
    basename = []
    snp_id = []
    snp_protein = []
    snps_db = []
    ensemble_protein_list = []
    bic_snps = []

    oms_in_vardb = []
    print "CSV\tDB\tOMS_in_BIC\tOMS_input\tOMS_out\tOMS_unique"
    for oms_in, bic_in in zip(oms_files, bic_files):
        db_lists = []
        #datastrucutres to store the mutations that are confirmed
        validated_identifier = []
        validated_sequence = []
        validated_tag = []
        validated_identifier_ens = []
        accepted_oms = []
        mass_diff = []
        temp_snps = []
        basename_oms = os.path.basename(oms_in)
        basename_bic = os.path.basename(bic_in)
        seq_dic, bic_mut_dic, bic_mut_tag_dic = read_biceps_log(bic_in)
        list_peps = []
        c = 0
        print basename_oms,
        print basename_bic,
        oms_results = pd.read_csv(oms_in, sep="\t")
        oms_results["unmod_sequence"] = [oms.AASequence(seqi).toUnmodifiedString() for seqi in oms_results["sequence"]]
        oms_results["single_uniprot"] = [BICl.get_uniprot(acci).replace("_MUT", "") for acci in oms_results["accessions"]]

        for seqi, uniproti in zip(oms_results["unmod_sequence"], oms_results["accessions"]):
            uniprot_ids_all = [BICl.get_uniprot(i).replace("_MUT", "") for i in uniproti.split(";") if i.count("|") >=2]
            # test if openms id in biceps results
            if seqi in bic_mut_dic:
                intersection_oms = np.unique(np.intersect1d(uniprot_ids_all, [uniprot_ids_seqi for uniprot_ids_seqi in bic_mut_dic[seqi].keys()]))
                if len(intersection_oms) >= 1:
                    c += 1
                    accepted_oms.append(True)

                    for item in intersection_oms:
                        tag = create_mutation_tag(item, seq_dic[seqi], seqi, uniprot_dic)
                        validated_identifier.append(item)
                        validated_sequence.append(seqi)
                        validated_tag.append(tag)
                        try:
                            validated_identifier_ens.append(mapping_dictionary[item])
                        except:
                            validated_identifier_ens.append("N.A")
                    list_peps.append(seq_dic[seqi].keys())
                    mass_diff.append(BICp.get_mass_diff(seqi, seq_dic[seqi].keys()[0]))
                else:
                    accepted_oms.append(False)
                    list_peps.append("")
                    mass_diff.append(0)
            else:
                    accepted_oms.append(False)
                    list_peps.append("")
                    mass_diff.append(0)


            # test if openms variation DB
            for uniprot_idi in uniprot_ids_all:
                if uniprot_idi in bic_mut_tag_dic:
                    pass
                else:
                    #probably a mutation with 2 residues
                    continue
                if uniprot_idi in mapping_dictionary:
                    ensemble_ids_temp =  mapping_dictionary[uniprot_idi][0].split(",")
                    for e_id in ensemble_ids_temp:
                        try:
                            # get more than ONE snp and transform to set
                            gt_snps = canprovar_df.loc[e_id]["native_id"].values
                        except:
                            # get more ONE snp and transform to set
                            gt_snps = [canprovar_df.loc[e_id]["native_id"]]
                            if len(gt_snps) != 1:
                                sys.exit("Wrong with the number of snps...")
                        #compare to BICEPS results
                        intersection = np.intersect1d(adjust_mutation_position(bic_mut_tag_dic[uniprot_idi]), gt_snps)

                        if len(intersection) != 0:
                            db_lists.append((e_id, ";".join(intersection), uniprot_idi))
                            temp_snps.extend(np.hstack(intersection))
                            snp_id.extend(intersection)
                            snp_protein.append(uniprot_idi)
                            ensemble_protein_list.append(e_id)
                            c += 1
        print c,
        snps_db.append(np.unique(temp_snps))
        fdr_c.append(c)
        basename.append(basename_oms)
        bic_snps.append(len(bic_mut_dic))

        oms_results["isBIC"] = accepted_oms
        oms_results["orig_seq"] = list_peps
        oms_results["mass_diff"] = mass_diff
        print oms_results.shape[0],
        oms_filtered = oms_results[(oms_results["isBIC"] == True) & (oms_results["mass_diff"] >= min_diff)]
        print oms_filtered.shape[0],
        oms_filtered.to_csv(out_ax+basename_oms+"_peptide_identifications.csv")
        print len(np.unique(oms_filtered["unmod_sequence"]))



        #create Provean data
        validated_df = pd.DataFrame()
        validated_df["identifier_uni"] = validated_identifier
        validated_df["identifier_ens"] = [i[0] if len(i) == 1 else i for i in validated_identifier_ens]
        validated_df["sequence"] = validated_sequence
        validated_df["tag"] = validated_tag
        validated_df.sort("identifier_uni", inplace=True)
        validated_df.drop_duplicates(inplace=True)
        write_to_provean_format(validated_df, out_ax+basename_oms+"_provean.csv", mass_filter=min_diff)

        #database variants
        df_vaars = pd.DataFrame()
        df_vaars["ensembl_id"] = [i[0] for i in db_lists]
        df_vaars["mutation_tags"] = [i[1] for i in db_lists]
        df_vaars["uniprot"] = [i[2] for i in db_lists]
        df_vaars = df_vaars.drop_duplicates()
        df_vaars.to_csv(out_ax+basename_oms+"_SNPs_canprovar.csv", sep="\t")

    #summary data frame for all files
    res_df = pd.DataFrame()
    res_df["SNPs_OMS"] = fdr_c
    res_df["SNPs_DB"] = [len(i) for i in snps_db]
    res_df["SNPs_BIC"] = bic_snps
    res_df["OMS_base"] = basename
    res_df["OMS"] = oms_files
    res_df["BIC"] = bic_files
    res_df.to_csv(out_ax+"SNPs_overview.csv", sep="\t")



if __name__ == "__main__":
   main()