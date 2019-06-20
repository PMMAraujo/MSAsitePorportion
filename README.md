# MSAsitePorportion

Python script to get the portotions of nucleotides or amino acids in each site of a multiple sequence alignment.


Usage:
  get_proportion_bilogical_sequences_by_site.py -i "input_file_name" -o "prefix_for-output_file_name" -c "conservancy_lvl"

Help:
  "-i","--input", required = True, help="Name of the input file, in fasta format only."
  "-o","--output", default = "output", type = str, help="Prefix for the name of the output files."
  "-n", "--nucleotides", action="store_true", help="Use when input is a nucleotides MSA - defaut is in aminoacids."
  "-c", "--conservancy_lvl", type = float, help="Lvl of conservancy desired for the outputs filtering."
