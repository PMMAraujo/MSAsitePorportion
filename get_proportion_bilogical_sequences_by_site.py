# -*- coding: utf-8 -*-
"""
Created on Tue Oct  3 20:13:16 2017

@author: Araujo
"""


import argparse
from Bio import AlignIO
import pandas as pd
import numpy as np

# Open file without prompt for test purposes
# alli = AlignIO.read("aas_teste.fasta", "fasta")

# Call user input of the desired multiple sequence alignement
#get_input = input("Please input sequence name: ")

# Pass input file to a biopython container for multiple sequence alignment
#alli = AlignIO.read(get_input, "fasta")


def Position_nucs(mutiple_alli, i):
    """ Get the nucleotide in a single position from all the sequences in the
    multiple sequence alignment."""
    nucleotides = []
    for record in mutiple_alli:
        nucleotides.extend(str(record.seq)[i].upper())
    return nucleotides

def Count_invalid(nucs_in_position):
    """ Get the number of characters that are nucleotides, not ambiguos."""
    not_nuc = set(nucs_in_position) - set(['A', 'C', 'D', 'E', 'F', 'G', 'H',
                                           'I', 'K', 'L', 'M', 'N', 'P', 'Q',
                                           'R', 'S', 'T', 'V', 'W', 'Y'])
    n_not_nuc = 0
    for i in not_nuc:
        n_not_nuc += nucs_in_position.count(i)
    return n_not_nuc


def Count_nuc(nucs_in_position):
    """ Count each nucleotide prevelance in a single position of the multiple
    sequence alignment. Using "Position_nucs" and "Count_invalid" functions
    outputs."""
    nucs_list = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N',
                 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    count_nucs = []
    length = len(nucs_in_position) - int(Count_invalid(nucs_in_position))
    for nuc in nucs_list:
        number = nucs_in_position.count(nuc)
        if int(length) == 0:
            percentage = 0
            count_nucs.append(percentage)
        else:
            percentage = int(number)/int(length)
            count_nucs.append(percentage)
    return count_nucs, length


def Nucleotides_per_site(mutiple_alli):
    """Creates a pandas dataframe (tsv file) being each raw a position of the
    multiple sequence alignment. Using "Count_nuc" function output"""
    data = pd.DataFrame([])
    length_max = len(mutiple_alli[0].seq)
    index = [str(x + 1) for x in range(int(length_max))]
    for i in range(length_max):
        counts = Count_nuc(Position_nucs(mutiple_alli, i))[0]
        real_length = Count_nuc(Position_nucs(mutiple_alli, i))[1]
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
#    df = data.sort_index(axis=1)
    return data


# Get resuling data frame into a avariable
#Result_df = Nucleotides_per_site(alli)


# Print data frame to file
#Result_df.to_csv('nucleotides_porportion_per_site.tsv', sep='\t')


#print("Job Done")

####
# Call user input of the desired multiple sequence alignement
#get_input = input("Please input the matrix: ")

# Pass input matrix to pnda dataframe and exclude "Valid nucleotides" column
#matrix = pd.DataFrame.from_csv(get_input, sep='\t')
#matrix = matrix.drop(['Valid nucleotides'], axis=1)

# Define threshold
#threshold = float(input("Please difine threshould value: "))

def Above_threshould(matrix, threshold):
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

#results = pd.DataFrame(np.array(Above_threshould(matrix, threshold)))
#results.columns = ['Position', 'Chracter', 'Value']
#results.to_csv("out_sum.csv", index=False, sep='\t')






if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--input", required = True,
                        help="name of the input file in fasta")
    parser.add_argument("-o","--output", default = "output", type = str,
                        help="name of the input file in fasta")    
    parser.add_argument("-c", "--conservancy_lvl", type = float,
                        help="lvl of conservancy desired for the outputs")

    args = parser.parse_args()
    
    if args.conservancy_lvl != None:
        print("")
        print("Building characters pooportion matrix and outputing the " 
        "postions above {0} conservancy "
        "treshould.".format(args.conservancy_lvl))        
        alli = AlignIO.read(args.input, "fasta")
        matrix = Nucleotides_per_site(alli)
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
              "-Positions above threshould: {1}".format(name_out_matrix, name_out_list))
        print("")
        print("Job Done")
    elif args.conservancy_lvl == None:
        print("")
        print("Only Building characters pooportion matrix.")
        alli = AlignIO.read(args.input, "fasta")
        Result_df = Nucleotides_per_site(alli)
        out_name = "matrix_{0}.tsv".format(args.output)
        Result_df.to_csv(out_name, sep='\t')
        print("")
        print("Outputs: \n -Porpotion Matrix: {0}".format(out_name))
        print("")
        print("Job Done")
    else:
        print("Error: invalid value for conservancy lvl")
