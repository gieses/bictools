# -*- coding: utf-8 -*-
"""
Created on Wed Oct 29 12:08:13 2014

@author: sven
"""
import re
import pandas as pd
import HTSeq
import matplotlib.pyplot as plt
import subprocess
import pyopenms as oms
import numpy as np
import brewer2mpl


class BICEPS_Hit():
    """
    Class to store BICEPS results. Explanation fo the single
    variables can be found on the github documentation:

    http://buotex.github.io/BICEPS/doc/html/index.html

    Paramters:
    --------------------------------------------------------
    spectrum, title, sequence, orig_sequence, fastaid, n, k, bic, penalty,
    score, tool, label, confidence)
    """

    def __init__(self, spectrum, title, sequence, orig_sequence, fastaid, n, k, bic, penalty, score, tool, label, confidence):
        self.spectrum = spectrum
        self.title = title
        self.sequence = sequence
        self.orig_sequence = orig_sequence
        self.fastaId = fastaid
        self.n = n
        self.bic = bic
        self.k = k
        self.penalty = penalty
        self.score = score
        self.tool = tool
        self.label = label
        self.confidence = confidence

    def fprint(self):
        """
        Adjusted print function.
        """
        print "spectrum:", self.spectrum
        print "title:", self.title
        print "sequence:", self.sequence
        print "orig_sequence:", self.orig_sequence
        print "fastaId:", self.fastaId
        print "n:", self.n
        print "bic:", self.bic
        print "k:", self.k
        print "penalty:", self.penalty
        print "score:", self.score
        print "tool:", self.tool
        print "label:", self.label
        print "confidence:", self.confidence


class BICEPS_Reader():
    """
    BICEPS Reader class.
    """
    def __init__(self, in_file):
        """
        Class to handle file objects from BICEPS results.

        Usage:
        ---------------------------------------
        bic_res = BICEPS_Reader(biceps_in)

        for PSM in bic_res:

            do something....
        """
        self.file = open(in_file, "r")

    def __iter__(self):
        # return each line transformed as list
        new_entry = True
        fasta_count = 1
        fasta_ids = []
        fasta_header = []
        for line in self.file:
            if new_entry is True:
                fasta_ids = []
                fasta_header = []
                spectrum = line.strip()
                new_entry = False
            if line.startswith("Title="):
                title = line[6:].strip()
            elif line.startswith("Sequence:"):
                sequence = re.search(': (\w+)', line).groups()[0]
            elif line.startswith("OrigSequence:"):
                orig_sequence = re.search(': (\w+)', line).groups()[0]
            elif line.startswith("fastaId:"):
                if fasta_count % 2 == 1:
                    fasta_ids.append(int(re.search(': (\d+)', line).groups()[0]))
                else:
                    fasta_header.append(line[7:].strip())
                fasta_count += 1
            elif line.startswith("n:"):
                n = int(re.search(': (\d+)', line).groups()[0])
            elif line.startswith("k:"):
                k = int(re.search(': (\d+)', line).groups()[0])
            elif line.startswith("bic:"):
                bic = float(line.split(":")[1])
            elif line.startswith("penalty:"):
                penalty = float(line.split(":")[1])
            elif line.startswith("score:"):
                score = float(line.split(":")[1])
            elif line.startswith("Tool:"):
                tool = line[6:].strip()
            elif line.startswith("Label:"):
                label = int(re.search(': (\d+)', line).groups()[0])
            elif line.startswith("Confidence:"):
                confidence = float(line.split(":")[1])
                new_entry = True
                yield BICEPS_Hit(spectrum, title, sequence, orig_sequence, zip(fasta_ids, fasta_header), n, k, bic, penalty, score, tool, label, confidence)


def get_uniprot(id_str):
    """
    Return the Uniprot identifier given a fasta header (whithout ">")

    Paramters:
    -------------------
    id_str: str,
            fasta sequence header
    """
    return(re.search(".*\|(.*)\|", id_str).groups()[0])


def read_uniprot_to_dic(in_file, mode="full"):
    """
    Reads a uniprot database and stores the identifier and sequence in a dic

    Paramters:
    ---------------------
        fasta_db: str,
                 file location for the fasta database

    mode: str,
          Either "full" or "seq". "seq" only stores the sequence while "full"
          stores the whole sequence object including name, description etc.

    Returns:
    -------------------------
    db_dic: dict,
            <key:value> with <uniprot_id>: Sequence
    """
    fasta = HTSeq.FastaReader(in_file)
    uniprot_dic = {}

    if mode == "seq":
        for seq in fasta:
            uniprot_dic[get_uniprot(seq.name)] = seq.seq
    elif mode == "full":
        for seq in fasta:
            uniprot_dic[get_uniprot(seq.name)] = seq
    else:
        print "Error! Unsupported mode: %s" % mode

    return(uniprot_dic)

    # -*- coding: utf-8 -*-
"""
Created on Wed Oct 29 14:59:09 2014

@author: sven
"""
import matplotlib.pyplot as plt
import subprocess
import re
import pyopenms as oms
import numpy as np
import pandas as pd
import brewer2mpl

def get_mass_diff(seq1, seq2):
    """
    Compute the mass difference for both sequences.
    """
    # full peptide, charge = 0
    m1 = oms.AASequence(seq1).getMonoWeight(0, 0)
    m2 = oms.AASequence(seq2).getMonoWeight(0, 0)
    return(np.abs(m1-m2))


def get_bic_cutoff(df, R_file, tmp_file, Rscipt_loc='/usr/bin/Rscript', fdr=0.05):
    """
    Retrieve the bic score cutoff.

    Paramters:
    ----------------------------------
    R_file : str,
             file location of the get_bic_cutoff file

    tmp_file : str,
             file location of a temporary file to use

    Rscipt_loc : str,
             Location of the the Rscript binary.
             E.g. /usr/bin/Rscript for linux. [default]

    Returns:
    ---------------------------------------
    bic_cutoff: float,
                bic score cut-off for FDR estimation at 5%

    """
    df.to_csv(tmp_file, sep="\t")
    # writes per default the results to
    # "C:\\tmp\\bic_test.csv"
    #'/home/sven/workspace/Spyder_workspace/py_coding/BICEPS/get_bic_cutoff.R'
    script_filename = R_file
    #'/home/sven/tmp/tmp.csv'
    param_filename = tmp_file+";{}".format(fdr)
    result = subprocess.Popen([Rscipt_loc, script_filename, param_filename], stdout=subprocess.PIPE)
    resis = result.communicate()[0]
    bic_cutoff = re.search("bic:(.*)", resis).groups()[0].strip()[0:-1]
    return(float(bic_cutoff))


def plot_bic_score(score_df, out_file, bic_score, file_type=".svg"):
    """
    Plots a histogram of the bic scores

    Paramters:
    --------------------------------
    score_df: df,
              df containing biceps results. Needs to have the two columns
              "group" and "bic"
    out_file: str,
              destination for the plotted figure
    file_type: str,
                output format [default: svg]

    bic_score: float,
                bic score cut-off for FDR estimation at 5%
    """
    #init plotting
    scores = []
    #group scores by yes/no containing subsitution
    group = score_df.groupby("group")

    #retrieve the cut-off
    for gr in group:
        scores.append(gr[1]["bic"].values)

    f, ax = plt.subplots(1, figsize=(11.69, 8.27))
    ax.hist(scores, label=["\wo subsitution", "\w subsitution"], bins=40)
    ax.set(xlabel="bic score", ylabel="Frequency")
    ax.grid(False)
    ax.axvline(bic_score, color="red", lw=3, label="%5 FDR cut-off ({:.2f})".format(bic_score))
    ax.legend(loc="upper left", ncol=1, prop={'size': 24})
    f.tight_layout()
    f.savefig(out_file+file_type, bbox_inches='tight', pad_inches=0.1)
    f.clf()


def plot_mass(in_file, out_file, sep="\t", file_type=".svg"):
    """
    Plot the mass difference between orig_sequence and sequence
    in the biceps results. in_file must have the two columns
    "sequence" and "orig_sequence".
    in_file = "/home/sven/data/BICEPS/HCT/results/HCT1/HCT1.csv"
    """
    df_bic = pd.read_csv(in_file, sep=sep)
    df_bic["mass_diff"] = [get_mass_diff(i.replace("Z","Q"), j.replace("Z","Q")) for i,j in zip(df_bic["orig_sequence"], df_bic["sequence"])]
    f, ax = plt.subplots(1, figsize=(11.69, 8.27))
    ax.boxplot(df_bic["mass_diff"])
    ax.set(xlabel="identifications", ylabel="mass difference", ylim=(0,100))
    ax.grid(False)
    ax.legend(loc="upper left", ncol=1, prop={'size': 24})
    f.tight_layout()
    f.savefig(out_file+file_type, bbox_inches='tight', pad_inches=0.1)
    f.clf()


def plot_mutation_matrix(in_file, out, file_type=".svg", norm=True):
    """
    Plots a heatmap for the mutations pattern. The input file is required
    to have the two columns "FROM" and "TO". Then these column combinations
    will be counted and plotted.

    Paramters:
    --------------------------
    in_file: str,
             dataframe source file
    out: str,
         outputfile
    file_type: str,
               file_type [ðefault:".svg"]
    norm: bool,
          True - normalize by all counts, False - don't normalize, show freqs

    in_file = "/home/sven/data/BICEPS/HCT/results/res_eval/provean.csv"
    """
    mut_dat = pd.read_csv(in_file, sep=" ")

    # init datastructures
    AA = ["A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V"]
    idx = np.arange(0, 20, 1)
    AA_map = {}
    #init the mapping dictionary
    for i,j in zip(AA, idx):
        AA_map[i] = j
    AA_mat = np.zeros(400).reshape(20,20)

    #iterate over dataframe and count mutations
    for FROM, TO in zip(mut_dat["FROM"], mut_dat["TO"]):
        row_i = AA_map[FROM]
        col_i = AA_map[TO]
        AA_mat[row_i][col_i] += 1

    # test
#    for aa1 in AA:
#        for aa2 in AA:
#            row_i = AA_map[aa1]
#            col_i = AA_map[aa2]
#            AA_mat[row_i][col_i] = get_mass_diff(aa1, aa2)
#    test_df = pd.DataFrame(AA_mat)
#    test_df.to_csv("/home/sven/data/BICEPS/HCT/test.csv", sep="\t")
    #normalize dataframe by sum?
    if norm:
        AA_mat = AA_mat / AA_mat.sum()
    else:
        pass
    #plot the figure
    #sns.set(context="poster",style="ticks", rc={"xtick.major.size": 8,"axes.labelsize": 14, "xtick.labelsize":14, "ytick.labelsize":14, "ytick.major.size": 8, "legend.fancybox":True})
    cmap = brewer2mpl.get_map("OrRd", "Sequential", 9).mpl_colormap
    plt.subplots(1, figsize=(11,8))
    plt.pcolormesh(AA_mat, cmap=cmap)
    plt.title("Preferred mutation pattern")
    plt.xlabel("FROM mutation")
    plt.ylabel("TO mutation")
    plt.yticks(idx+0.5, AA)
    plt.xticks(idx+0.5, AA)
    plt.colorbar()
    plt.savefig(out+"pref_mutation"+file_type)
    plt.tight_layout()
    plt.savefig(out+file_type, bbox_inches='tight', pad_inches=0.1)
    plt.clf()