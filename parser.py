# -*- coding: utf-8 -*-
"""
Created on Thu Oct  5 15:09:03 2017

@author: Araujo
"""
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--input", required = True,
                        help="name of the input file in fasta")
    parser.add_argument("-o","--output", default = "output", type = str,
                        help="name of the input file in fasta")
    parser.add_argument("-n", "--nucleotides", action="store_true",
                        help="input in nucleotides - defaut is in aminoacids")
    parser.add_argument("-c", "--conservancy_lvl", type = float,
                        help="lvl of conservancy desired for the outputs")

    args = parser.parse_args()

    if args.nucleotides == True:
        print("nucleotides")
    else:
        print("aminoacids")
    
    
    
    
"""    
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
              "-Positions above threshould: {1}".format(name_out_matrix,
                                                        name_out_list))
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
        
"""