# -*- coding: utf-8 -*-
"""
Created on Wed Oct 29 11:59:16 2014

@author: sven
"""

import argparse
import textwrap
import pandas as pd
import bictools as BIC


def bic_to_list(in_file):
    """
    Function that transforms the BICEPS results to a dataframe. The dataframe
    will have all columns present in the BICEPS result file.
    See http://buotex.github.io/BICEPS/doc/html/index.html
    """

    #init file object
    bic_res = BIC.BICEPS_Reader(in_file)

    # store results
    PSM_list = []
    for PSM in bic_res:
        PSM_list.append(PSM)
    return(PSM_list)


def to_dataframe(bic_results):
    """
    Transform list of biceps results to a dataframe. Add a column that
    states if the PSM contains a SNP or not.
    """
    df = pd.DataFrame()
    df["spectrum"] = [i.spectrum for i in bic_results]
    df["title"] = [i.title for i in bic_results]
    df["sequence"] = [i.sequence for i in bic_results]
    df["orig_sequence"] = [i.orig_sequence for i in bic_results]
    df["n"] = [i.n for i in bic_results]
    df["bic"] = [i.bic for i in bic_results]
    df["k"] = [i.k for i in bic_results]
    df["penalty"] = [i.penalty for i in bic_results]
    df["score"] = [i.score for i in bic_results]
    df["tool"] = [i.tool for i in bic_results]
    df["label"] = [i.label for i in bic_results]
    df["confidence"] = [i.confidence for i in bic_results]
    df["fastaId"] = [i.fastaId for i in bic_results]
    df["isSNP"] = [1 if i.sequence.upper() != i.orig_sequence.upper() else 0 for i in bic_results]

    return(df)

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

    parser.add_argument('--infile', metavar='infile', type=str, nargs='+',
                       help='BICEPS input file', required=True)
    parser.add_argument('--outfile', metavar='outfile', type=str, nargs='+',
                       help='Desired outputfile in .csv format', required=True)

    args = parser.parse_args()
    infile = args.infile[0]
    outfile = args.outfile[0]

    print ""
    print "###################################################################"
    print "# input: %s" % infile
    print "#"
    print "# output: %s" % outfile
    print "###################################################################"
    print ""

    #transform to peptide list
    PSM_list = bic_to_list(infile)

    #transform to dataframe
    df = to_dataframe(PSM_list)
    print "%s PSMs found in input file." % len(PSM_list)
    print "%s PSMs contain sequence variations." % df["isSNP"].sum()
    df.to_csv(outfile, sep="\t")


#entry point
if __name__ == "__main__":
    main()
