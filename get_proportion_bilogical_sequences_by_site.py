# -*- coding: utf-8 -*-
"""
Created on Tue Oct  3 20:13:16 2017

@author: Araujo
"""

import argparse
from Bio import AlignIO
import pandas as pd
import numpy as np


def This_pos_chars(mutiple_alli, i):
    """ Get the characters (nucleotides or aminoacids) in a single position 
    from all the sequences in the multiple sequence alignment."""
    nucleotides = []
    for record in mutiple_alli:
        nucleotides.extend(str(record.seq)[i].upper())
    return nucleotides
            
            
def Count_invalid(nucs_in_position):
    """ Count number of characters that are invalid: gaps, ambiguities or 
    strange characters."""
    not_nuc = set(nucs_in_position) - set(characters_list)
    n_not_nuc = 0
    for i in not_nuc:
        n_not_nuc += nucs_in_position.count(i)
    return n_not_nuc


def Count_char(nucs_in_position):
    """ Count each character (20 structural aas or 4 nucleotides) prevelance 
    in a single position of the multiple sequence alignment. Using 
    "This_pos_chars" and "Count_invalid" functions outputs."""
    count_nucs = []
    length = len(nucs_in_position) - int(Count_invalid(nucs_in_position))
    for nuc in characters_list:
        number = nucs_in_position.count(nuc)
        if int(length) == 0:
            percentage = 0
            count_nucs.append(percentage)
        else:
            percentage = int(number)/int(length)
            count_nucs.append(percentage)
    return count_nucs, length

                                         
def Nucleotides_per_site(mutiple_alli):
    """Creates a pandas dataframe (tsv file) being each row a position of the
    multiple sequence alignment, and each column a nucleotide.
    Using "Count_char" function output"""
    data = pd.DataFrame([])
    length_max = len(mutiple_alli[0].seq)
    index = [str(x + 1) for x in range(int(length_max))]
    for i in range(length_max):
        counts = Count_char(This_pos_chars(mutiple_alli, i))[0]
        real_length = Count_char(This_pos_chars(mutiple_alli, i))[1]
        data = data.append(pd.DataFrame({'A': [counts[0]], 'C': [counts[1]],
                                         'G': [counts[2]], 'T': [counts[3]],
                                         'Valid nucleotides': [real_length]}))
    data.index = index
    data.index.name = "Position"
    return data

    
def Aas_per_site(mutiple_alli):
    """Creates a pandas dataframe (tsv file) being each row a position of the
    multiple sequence alignment, and each column an aminoacid. Using 
    "Count_char" function output"""
    data = pd.DataFrame([])
    length_max = len(mutiple_alli[0].seq)
    index = [str(x + 1) for x in range(int(length_max))]
    for i in range(length_max):
        counts = Count_char(This_pos_chars(mutiple_alli, i))[0]
        real_length = Count_char(This_pos_chars(mutiple_alli, i))[1]
        data = data.append(pd.DataFrame({'A': [counts[0]], 'C': [counts[1]],
                                         'D': [counts[2]], 'E': [counts[3]],
                                         'F': [counts[4]], 'G': [counts[5]],
                                         'H': [counts[6]], 'I': [counts[7]],
                                         'K': [counts[8]], 'L': [counts[9]],
                                         'M': [counts[10]], 'N': [counts[11]],
                                         'P': [counts[12]], 'Q': [counts[13]],
                                         'R': [counts[14]], 'S': [counts[15]],
                                         'T': [counts[16]], 'V': [counts[17]],
                                         'W': [counts[18]], 'Y': [counts[19]],
                                         'Valid nucleotides': [real_length]}))
    data.index = index
    data.index.name = "Position"
    return data


def Above_threshould(matrix, threshold):
    """ In case of a threshould being specified ( option "-c") this function
    creates a list of all the characters with proportion above threshould
    for its alignment position. Starts from "Aas_per_site" or 
    "Nucleotides_per_site" function outputs."""
    matrix = matrix.drop(['Valid nucleotides'], axis=1)
    results = []        
    for index, row in matrix.iterrows():
        for column in matrix:
            value = matrix.loc[index,matrix[column].name]
            if value >= threshold:
                results.append([index, matrix[column].name, value])
            else:
                pass
    return results

    
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--input", required = True,
                        help="Name of the input file, in fasta format only.")
    parser.add_argument("-o","--output", default = "output", type = str,
                        help="Prefix for the name of the output files.")
    parser.add_argument("-n", "--nucleotides", action="store_true",
                        help="Use when input is a nucleotides MSA - defaut "
                        "is in aminoacids.")
    parser.add_argument("-c", "--conservancy_lvl", type = float,
                        help="Lvl of conservancy desired for the outputs "
                        "filtering.")
    args = parser.parse_args()
    
# Based on datatype (option -n) decide which function to use to build the matrix
    if args.nucleotides == True:
        print("")
        print("Data type: Nucleotides")
        alli = AlignIO.read(args.input, "fasta")
        print("Building characters poportion matrix.")
        characters_list = ['A', 'C', 'G', 'T']
        matrix = Nucleotides_per_site(alli)
    elif args.nucleotides == False:
        print("")
        print("Data type: Aminoacids")   
        alli = AlignIO.read(args.input, "fasta")
        print("Building characters pooportion matrix.")
        characters_list = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
                           'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
        matrix = Aas_per_site(alli)
    else:
        print("Error: could not understand data type option (-n).")
      
# If threshould if difined a list above treshould is outputed with the matrix       
    if args.conservancy_lvl != None:
        print("")
        print("Outputing the postions above {0} conservancy "
        "treshould.".format(args.conservancy_lvl))        
        threshold = args.conservancy_lvl
        above = Above_threshould(matrix, threshold)
        results = pd.DataFrame(np.array(above))
        results.columns = ['Position', 'Chracter', 'Value']
        
        name_out_matrix = "matrix_{0}.tsv".format(args.output)
        name_out_list = "above_threshould_{0}.tsv".format(args.output)
        matrix.to_csv(name_out_matrix, sep='\t')
        results.to_csv(name_out_list, index=False, sep='\t')
        print("")
        print("Outputs: \n -Porpotion Matrix: {0}\n "
              "-Positions above threshould: {1}".format(name_out_matrix,
                                                        name_out_list))
        print("")
        print("Job Done")
    elif args.conservancy_lvl == None:
        out_name = "matrix_{0}.tsv".format(args.output)
        matrix.to_csv(out_name, sep='\t')
        print("")
        print("Outputs: \n -Porpotion Matrix: {0}".format(out_name))
        print("")
        print("Job Done")
    else:
        print("Error: invalid value for conservancy lvl")
        