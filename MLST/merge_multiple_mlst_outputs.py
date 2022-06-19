import subprocess as sp
import os
import threading
from types import new_class
import numpy as np
import sys


# folder containing only the mlst files as command line argument
# /media/milena/One\ Touch/mlst/jejuni_mlst_all_outfiles

mlst_directory = sys.argv[1] 

code_dir = os.getcwd()
os.chdir(mlst_directory)
print('\nDirectory with the mlst out files:')
print(os.getcwd())
all_mlst_files = os.listdir()
#only files with .out ending
all_mlst_files = [i for i in all_mlst_files if '.out' in i]

# parse one line in the mlst form the mlst out format to the format supported by philovyz
# The new_st_type argument is the ST assigned to a genome when the mlst program didn't assign one on its own
# It has to be counted outside the function
def parse_genome_line(line, new_st_count):
	accession, schema, st, gene1, gene2, gene3, gene4, gene5, gene6, gene7 = line.split('\t')
	# remove filepath from accession numbers
	accession = "GCA_"+ accession.split("GCA_")[1]
	accession = accession[0:13]
	# extract allele types
	genes = [gene1, gene2, gene3, gene4, gene5, gene6, gene7]
	allele_numbers = [i.split('(',1)[1] for i in genes] #remove front (
	allele_numbers = [i.split(')',1)[0] for i in allele_numbers] #remove back )
	# merge allele types and sequence numbers
	# add sequence type if none is assigned
	if st == "-": #when the st is unknown by the mlst program
		all_numbers = str(new_st_count)+"\t"+"\t".join(allele_numbers)
		new_st_count +=1
	else:
		all_numbers = st+"\t"+"\t".join(allele_numbers)
	return [accession, all_numbers]

# open merged output file:
with open('mlst_all_ST.tsv', 'w') as st_out, open('metadata.txt', 'w') as meta_out, open('combined_mlst.txt', 'w') as combined_out:
	
	# open first file to get the header
	print('file used for header: '+all_mlst_files[0])
	header_file = open(all_mlst_files[0], 'r')
	lines = header_file.readlines()
	accession, schema, st, gene1, gene2, gene3, gene4, gene5, gene6, gene7 = lines[0].split('\t')
	genes = [gene1, gene2, gene3, gene4, gene5, gene6, gene7]
	gene_names = [i.split('(',1)[0] for i in genes] # leave out numbers from first line
	header = "ST\t" + "\t".join(gene_names)
	st_out.write(header+"\n")
	combined_out.write("Accession_nr\t"+header+"\n")
	meta_out.write("Accession_nr\n")
	header_file.close()
	
	#counting index for new sequence types (high enough to not interfere with preexisting STs)
	new_st_count = 50000 
	
	
	
	for mlst_file in all_mlst_files:
		print(mlst_file)
		
		with open(mlst_file, 'r') as mlst_in:
			lines = mlst_in.readlines()
			print('no. of genomes in file: ', len(lines[1:]))
			for line in lines[1:]:
				accession, all_numbers = parse_genome_line(line, new_st_count)
				new_st_count +=1
				if "," in all_numbers or "~" in all_numbers:
					pass
				else: 
					st_out.write(all_numbers+"\n")
					combined_out.write(accession+"\t"+all_numbers+"\n")
					meta_out.write(accession+"\n")
	
os.chdir(code_dir)
					
