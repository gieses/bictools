# -*- coding: utf-8 -*-
"""
Created on Wed Oct 29 14:19:58 2014

@author: sven
"""


import argparse
import textwrap
import pandas as pd
import bictools as BICl
import bictools as BICp
import re
import glob
import subprocess
import HTSeq
import os
import numpy as np
import matplotlib.pyplot as plt
import sys

#==============================================================================
# WRITE FUNCTIONS
#==============================================================================

def write_log(out_file, collect_mutations):
    """
    Write a log file for the sequences that were found in BICEPS with
    substitutions.

    Parameters:
    -----------------------------
    out_file: str,
                outfile anme

    collect_mutations: dic,
                       dictionary containing data

    Format:
    -------------
    'Q9Y6J6': [('T9S;E11D', 'MSTLSNFTQTLEDVFR', 'MSTLSNFTQSLDDVFR')],
    'Q9Y6R4': [('A754S', 'QAGKLFCDIAGMLLK', 'QAGKLFCDISGMLLK')],
    """
    print "write log...",
    mut_log = open(out_file, "w")
    mut_log.write("uniprot_id,tag,orig_sequence,sequence\r\n")
    counter = 0
    counter_2 = 0
    for uni_id, mutation_data in collect_mutations.iteritems():
        tag_dic = {}
        for i, tagi in enumerate(mutation_data):
            if tagi[0] in tag_dic:
                counter_2 += 1
                pass
            else:
                tag_dic[tagi[0]] = ""
                #print("%s,%s,%s,%s\r\n" % (uni_id, tagi[0], tagi[1], tagi[2]))
                mut_log.write("%s,%s,%s,%s\r\n" % (uni_id, tagi[0], tagi[1], tagi[2]))
                counter += 1
    mut_log.close()


def write_whitelist(out_file, collect_mutations):
    """
    Write a log file for the sequences that were found in BICEPS with
    substitutions.

    Parameters:
    -----------------------------
    out_file: str,
                outfile anme

    collect_mutations: dic,
                       dictionary containing data

    Format:
    -------------
    'Q9Y6J6': [('T9S;E11D', 'MSTLSNFTQTLEDVFR', 'MSTLSNFTQSLDDVFR')],
    'Q9Y6R4': [('A754S', 'QAGKLFCDIAGMLLK', 'QAGKLFCDISGMLLK')],
    """
    print "write log...",
    #mut_log = open(path+"/"+basename+".mut", "w")
    whitelist = open(out_file, "w")
    counter = 0
    counter_2 = 0
    for uni_id, mutation_data in collect_mutations.iteritems():
        #dictionary to keep track of tacsk already written
        tag_dic = {}
        for i, tagi in enumerate(mutation_data):
            if tagi[0] in tag_dic:
                counter_2 += 1
                pass
            else:
                tag_dic[tagi[0]] = ""
                counter += 1
                whitelist.write(">sp|%s_MUT %s\r\n" % (uni_id, tagi[0]))
                whitelist.write("%s\r\n" % (tagi[2]))

    whitelist.close()


def write_db_rd(out_file, collect_mutations, uniprot_dic, append_UP_seqs=False):
    """
    Write a new FASTA target database file for the sequences that were found in
    BICEPS with substitutions. Writes for each SNP a new entry.
    Only the BICEPS sequences are added to the target database. Other proteins
    will not be added.

    Parameters:
    -----------------------------
    out_file: str,
                outfile anme

    collect_mutations: dic,
                       dictionary containing data

    Format:
    -------------
    'Q9Y6J6': [('T9S;E11D', 'MSTLSNFTQTLEDVFR', 'MSTLSNFTQSLDDVFR')],
    'Q9Y6R4': [('A754S', 'QAGKLFCDIAGMLLK', 'QAGKLFCDISGMLLK')],

    uni_id = "P30613"
    mutation_data = collect_mutations[uni_id]
    """
    print "write db...",
    #mut_log = open(path+"/"+basename+".mut", "w")
    whitelist = open(out_file, "w")
    counter = 0
    counter_2 = 0
    test_db = {}
    db_logger = {}
    #test_dic = {}
    for uni_id, mutation_data in collect_mutations.iteritems():
        #dictionary to keep track of tacsk already written
        tag_dic = {}
        db_logger[uni_id] = ""
        if uni_id == "P13929":
            pass
        for i, tagi in enumerate(mutation_data):
            if tagi[0] in tag_dic:
                counter_2 += 1
                pass
            else:
                tag_dic[tagi[0]] = ""
                counter += 1

        #for each mutation create a new fasta entry
        for mutation_id in tag_dic:

            temp_sequence = list(uniprot_dic[uni_id])
            # if two mutations are found in one spectrum iterate through
            # both of them
            mtags = mutation_id.split(";")
            for mtagsi in mtags:
                FROM, TO, POS = split_mutation_id(mtagsi)
                temp_sequence[POS] = TO
            whitelist.write(">sp|%s_MUT|%s %s\r\n" % (uni_id, mutation_id, mutation_id))
            whitelist.write("%s\r\n" % ("".join(temp_sequence)))
            #print(">sp|%s_MUT|%s %s\r\n" % (uni_id, mutation_id, mutation_id))
            #print("%s\r\n" % ("".join(temp_sequence)))

    #append normal uniprot sequences to the database
    # TODO:
    if append_UP_seqs:
        for uni_key in uniprot_dic.keys():
            if uni_key in db_logger:
                pass
            else:
                whitelist.write(">sp|%s| \r\n" % (uni_id, mutation_id, mutation_id))
                whitelist.write("%s\r\n" % uniprot_dic[uni_key])

    print "(%s/%s redundant tags)" % (counter_2, counter),
    whitelist.close()
    return(test_db)


#==============================================================================
# HELPER FUNCTIONS
#==============================================================================

def split_mutation_id(mutation_tag):
    """
    Splits the mutation tag into FROM, TO, POS.

    Paramters:
    --------------------
    mutation_tag: str,
                    mutation tag, e.g. "D192A"

    Returns:
    -------------------
    FROM, TO, POS: char, char, int,
                    e.g. D, A, 192
    """

    FROM, POS, TO = re.search("([A-Z*-]+)(\d+)([A-Z*-]+)", mutation_tag).groups()
    POS = int(POS)
    return (FROM, TO, POS)


def get_difference(original_sequence, sequence):
    """
    Returns the FROM to Mutation for the given sequences.

    Parameters:
    ----------------------------
    original_sequence: str,
             AASequence from Biceps result as found in the initial DB
    sequence: str,
             AASequence from Biceps results as identified by the MGF
    """
    #we can have more than two subsitutions
    FROM = []
    TO = []
    POS = []
    for i, item in enumerate(original_sequence):
        if original_sequence[i].upper() != sequence[i].upper():
            FROM.append(original_sequence[i])
            TO.append(sequence[i])
            POS.append(i)
    return(FROM, TO, POS)


def create_mutation_tag(FROM, TO, POS, original_peptide, protein_sequence):
    """
    Function to create mutation tags. The function requires the
    identified peptide, the original peptide sequence as well as the
    protein sequence.

    Paremter:
    ----------------------------
    FROM: char,
          character that was initial found in the DB
    TO: char,
         character (AA) was was used to replace the initial DB entry
    POS: int,
        the posiiton in the sequence
    original_peptide: str,
                amino acid peptide
    protein_sequence: str,
                proteins equence in AA

    Returns:
    ------------------------
    tags: list,
             list of uniprot identifier with mutation tag to run the online
             mapping service ()
    FROM, TO, POS, original_peptide, protein_sequence = \
    FROM, TO, POS, PH.orig_sequence, fasta_DB[uniprot_id].seq
    """
    pep_index = protein_sequence.index(original_peptide)
    POS = np.array(POS)
    pos = pep_index + POS

    tags = []
    for fi, ti, posi in zip(FROM, TO, pos):
        tags.append("%s%s%s" % (fi, posi, ti))
        #print fi, ti, posi
        #print protein_sequence[posi]
    return(";".join(tags))


def compute_cut_off_single(biceps_in, R_file, tmp_file, rscript_file, fdr):
    """
    Computes the BIC score cut-off to disgard certain subsitutions before
    the creation of the variant database

    Paramter:
    -------------------------
    biceps: str,
            text file origin


    Returns:
    -------------------------
    cut_off: float,
             5% FDR cut_off
    """

    bic_res = BICl.BICEPS_Reader(biceps_in)
    df = pd.DataFrame()
    #iterate over biceps results and only deal with cases
    #where the sequences are NOT equal
    bic = []
    for PH in bic_res:
        bic.append(PH.bic)

    df["bic"] = bic
    cut_off = BICp.get_bic_cutoff(df, R_file, tmp_file, rscript_file, fdr=fdr)
    return(cut_off)

def main():
    """
    Main function
    """

    parser = argparse.ArgumentParser(
         prog='bic_DB.py',
         formatter_class=argparse.RawDescriptionHelpFormatter,
         description=textwrap.dedent('''\

         ---------------------------------------------------------------------
         | Description:                                                      |
         ---------------------------------------------------------------------
         Command line tool to create customized databases from
         biceps results. ("\\t").

         The input requires the biceps search results, the target databases
         and a desired output directory. The biceps search results can either
         be a directory or a single file.


         ---------------------------------------------------------------------
         | Usage:                                                            |
         ---------------------------------------------------------------------
         python bic_DB.py --in_bic Biceps.Results.HCT_all.txt \\
                          --in_db HS_reviewed_19092014.fasta \\
                          --out_dir /home/results


         =====================================================================
                 '''))

    parser.add_argument('--in_bic', metavar='in_bic', type=str,
                       help='BICEPS input file directory or file', required=True)
    parser.add_argument('--in_db', metavar='in_db', type=str,
                       help='protein database used for the --in_bic file', required=True)
    parser.add_argument('--out_dir', metavar='out_dir', type=str,
                       help='Desired outputfile in .fasta format', required=True)

    parser.add_argument('--rscript', metavar='rscript', type=str,
                       default="/usr/bin/Rscript",
                       help='Location of the Rscript binary.')
    parser.add_argument('--tmp_file', metavar='tmp_file', type=str,
                       default="/home/sven/tmp/tmp.csv",
                       help='A file that is used for writing temporary files')
    parser.add_argument('--R_file', metavar='R_file', type=str,
                       default="",
                       help='Location of the R script to retrieve the cut-off\
                       for confidence estimation. If unset the score\
                       distributions will not be plotted.')
    parser.add_argument('--mass_diff', metavar='mass_diff', type=float,
                       help='Mass filter for minimum difference between the \
                       sequence identified by biceps and the original sequence', default=-1.0
                       )
    parser.add_argument('--bic_fdr', metavar='bic_fdr', type=bool,
                        choices=[1, 0],
                       help='Use the bic FDR filter to limit the number \
                       of variant candidates in the first instance. Only \
                       peptides above the FDR cut-off will be used to create \
                       the variant database.', default=False
                       )
    parser.add_argument('--max_mut', metavar='max_mut', type=int, default=1,
                       help='Allowed number of mutations for peptides \
                       the be reported in the candidate database.')
    parser.add_argument('--fdr', metavar='fdr', type=float, default=0.05,
                       help='FDR threshold to estimate via CurveFDP.')
    #==============================================================================
    # Example
    #==============================================================================
    # /home/sven/workspace/Spyder_workspace/py_coding/BICEPS/get_bic_cutoff.R
#    in_bic = "/home/sven/data/BICEPS/HCT/results/Biceps.Results.HCT_all_1.txt"
#    in_db  = "/home/sven/data/BICEPS/HCT/HS_reviewed_19092014.fasta"
#    out_dir ="/home/sven/test/biceps/"
#    tmp_file = "/home/sven/tmp/"
#    min_mass_diff = -1
#    tmp_file = "/home/sven/tmp/tmp.csv"
#    R_file = "/home/sven/workspace/Spyder_workspace/py_coding/BICEPS/get_bic_cutoff.R"
#    rscript_file = "/usr/bin/Rscript"


    args = parser.parse_args()
    in_bic = args.in_bic
    in_db = args.in_db
    out_dir = args.out_dir
    min_mass_diff = args.mass_diff
    max_mut = args.max_mut
    fdr = args.fdr

    if args.R_file == "":
        do_plot = False
    else:
        do_plot = True
        rscript_file = args.rscript
        tmp_file = args.tmp_file
        R_file = args.R_file
        bic_fdr = args.bic_fdr

    print ""
    print "###################################################################"
    print "# db input: %s" % in_db
    print "# biceps input: %s" % in_bic
    print "#"
    print "# output dir: %s" % out_dir
    print "# min mass diff: %s" % min_mass_diff
    print "# max mutations: %s" % max_mut
    print "# bic_FDR enabled: %s" % args.bic_fdr
    print "###################################################################"
    print ""

    #==========================================================================
    #     START MAIN PROGRAM
    #==========================================================================
    #test values
    #in_db = '/home/sven/data/BICEPS/fasta/HS_reviewed_19092014.fasta'
    #out_dir = "/home/sven/data/BICEPS/HeLa/results/HCT1/out/"
    #in_bic = "/home/sven/data/BICEPS/HCT/results/HCT1/"

    #init list of input files
    if os.path.isfile(in_bic):
        print ("single file mode:")
        bic_files = [in_bic]
    elif os.path.isdir(in_bic):
        bic_files = sorted(glob.glob(in_bic+"Biceps.Results.*.txt"))
        print ("batch mode: %s files" % len(bic_files))
    else:
        print "Error! -in_bic is not a correct file or directory!"
        sys.exit()

    if not os.path.exists(out_dir):
        print "%s created!" % out_dir
        os.makedirs(out_dir)

    #==========================================================================
    #     start actual analysis
    #==========================================================================
    # id-sequence mapping
    # Example:
    # 'Q8TBF5' - 'MAARVAAVRAAAWLLLGAATGLTRGPAAAFTAAR....."
    uniprot_dic = BICl.read_uniprot_to_dic(in_db, mode="seq")
    # 'Q8TBF5' -  HTSeq object with sequence, name, description
    fasta_DB = BICl.read_uniprot_to_dic(in_db, mode="full")

    #some summary datastructures
    file_list = []
    hits_ws_list = []
    hits_wos_list = []
    ov_cutoff = []
    mass_diff_counter = 0
    print "analyse biceps results..."
    for file_counter, biceps_in in enumerate(bic_files):

        #first get the confidence cut-off if option is set
        if bic_fdr and do_plot:
            #reuse cut_off computation
            bic_cut_off = compute_cut_off_single(biceps_in, R_file, tmp_file, rscript_file, fdr)
            cut_off = bic_cut_off
        elif bic_fdr:
            #only for FDR cut_off
            bic_cut_off = compute_cut_off_single(biceps_in, R_file, tmp_file, rscript_file, fdr)
        elif do_plot:
            # only for plotting
            cut_off = compute_cut_off_single(biceps_in, R_file, tmp_file, rscript_file, fdr)
            bic_cut_off = -10000
        else:
            #set dummy cut-off to reuse code
            bic_cut_off = -10000

        print "%s/%s" % (file_counter, len(bic_files))
        scores = []
        label = []
        # dictionary for hits with changed amino acids
        mutable_DB = BICl.read_uniprot_to_dic(in_db, mode="full")

        #mutation dictionary
        bic_mut_dic = {}

        #counter for some stats
        hits_wos = 0
        hits_ws = 0

        biceps_mut_dic = {}
        modified_uniprots = {}
        # set current filename (without path) and path seperately
        basename = os.path.basename(biceps_in)


        #init whitelist for biceps sequences

        whitelist_db_file_seq = open(out_dir+"/whitelist_"+basename+"_seq.fasta", "w")

        #new fasta file
        #used as input db for nomral search algorithms!
        new_fasta_db_file = open(out_dir+"/SAAV_db_"+basename+".fasta", "w")

        #mutation dic
        collect_mutations = {}

        #init biceps iterator
        bic_res = BICl.BICEPS_Reader(biceps_in)

        lower_cut_off_counter_saav = 0
        upper_cut_off_counter_saav = 0

        #iterate over biceps results and only deal with cases
        #where the sequences are NOT equal
        for PH in bic_res:
            PH.sequence = PH.sequence.upper()
            PH.orig_sequence = PH.orig_sequence.upper()
            if PH.sequence == PH.orig_sequence:
                hits_wos += 1
                scores.append(PH.bic)
                label.append("0")
                continue
            else:
                #score log
                scores.append(PH.bic)
                label.append("1")

                #hit has enough confidence?
                if PH.bic >= bic_cut_off:
                    upper_cut_off_counter_saav += 1
                else:
                    lower_cut_off_counter_saav += 1
                    continue

                #hit has large enough mass difference?
                # e.g. K <--> Q mutations have a very low mass difference
                # one could disable the report of these mutations
                if BICp.get_mass_diff(PH.sequence, PH.orig_sequence) >= min_mass_diff:
                    #prepare mutation identifier
                    FROM, TO, POS = get_difference(PH.orig_sequence, PH.sequence)

                    #continue if number of mutations is too high
                    # double mutations are often cumbersome...
                    # reverted mutations or missing fragments
                    # allow for recombination of amino acids to achieve a given
                    if len(FROM) > max_mut:
                        continue

                    hits_ws += 1
                    uniprot_ids = [BICl.get_uniprot(i[1]) for i in PH.fastaId]
                    for uniprot_id in uniprot_ids:
                        tags = create_mutation_tag(FROM, TO, POS, PH.orig_sequence, fasta_DB[uniprot_id].seq)
                        if uniprot_id in collect_mutations:
                            collect_mutations[uniprot_id].extend([(tags, PH.orig_sequence, PH.sequence)])
                        else:
                            collect_mutations[uniprot_id] = [(tags, PH.orig_sequence, PH.sequence)]
                else:
                    mass_diff_counter += 1
        #write file with mutation tags, uniprot id, from_peptide, to_peptide
        write_log(out_dir+"/"+basename+".mut", collect_mutations)
        print "\t", out_dir+"/"+basename+".mut :"

        #write file with uniprot id + mutated peptide
        write_whitelist(out_dir+"/whitelist_"+basename+"_seq.fasta", collect_mutations)
        print "\t", out_dir+"/whitelist_"+basename+"_seq.fasta :"


        #bic_db = write_db(path+"/mod_db_"+basename+".fasta", collect_mutations, uniprot_dic)
        #test_bic_dic(bic_db, path+"/mod_db_"+basename+".fasta", collect_mutations, uniprot_dic)

        #write redundant database form
        write_db_rd(out_dir+"/SAAV_db_"+basename+".fasta", collect_mutations, uniprot_dic)
        #write the reference fasta database with the full protein sequence
        print "\t", out_dir+"/SAAV_db_"+basename+".fasta"
        score_df = pd.DataFrame()
        score_df["bic"] = scores
        score_df["group"] = label
        if do_plot:
            BICp.plot_bic_score(score_df, out_dir+basename+"_bic_scores", cut_off)
            ov_cutoff.append(cut_off)
        else:
            cut_off = np.nan

        hits_ws_list.append(hits_ws)
        hits_wos_list.append(hits_wos)
        print ""
        print "%s ..... DONE!:" % basename,
        whitelist_db_file_seq.close()
        new_fasta_db_file.close()

        print "Ignored {} variant candidates because of similar peptide \
        masses".format(mass_diff_counter)

        #write a short summary
        fobj = open(out_dir+basename+"_summary.txt", "w")
        fobj.write("""\r\n \
Program Call: {0} \r\n \
basefile: {1}\r\n \
bic cut_off: {2} \r\n \
hits \wo substitution: {3} \r\n \
hits \w substitution: {4} \r\n \
hits ignored (min mass difference): {5} \r\n \
SAAV candidates < bic cut_off: {6} \r\n \
SAAV candidates >= bic cut_off: {7} \r\n \
""".format(args, basename, cut_off, hits_wos,
                   hits_ws, mass_diff_counter, lower_cut_off_counter_saav,
                   upper_cut_off_counter_saav))
        fobj.close()


#entry point
if __name__ == "__main__":
    main()
