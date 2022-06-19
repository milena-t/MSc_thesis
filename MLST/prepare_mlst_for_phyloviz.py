import subprocess as sp
import os
import threading
from types import new_class
import numpy as np
import sys

# /home/milena/masters_thesis/thesis_git/mlst/jejuni1000.mlst
# edited mlst file as command line argument
# edit with find/replace to:
#  remove the gene names from all the cells so that there are only numbers left, 
#  remove the ~ 
# the arguments are organized like: python3 mlst.py[0] file[1]

code_dir = os.getcwd()

mlst_original = sys.argv[1] # [0] would be the name of this file



with open(mlst_original, 'r') as mlst_infile, open('sequence_types.txt', 'w') as st_out, open('metadata.txt', 'w') as meta_out, open('combined_mlst.txt', 'w') as combined_out:
    lines = mlst_infile.readlines()

    # make header for st_out
    accession, schema, st, gene1, gene2, gene3, gene4, gene5, gene6, gene7 = lines[0].split('\t')
    genes = [gene1, gene2, gene3, gene4, gene5, gene6, gene7]
    gene_names = [i.split('(',1)[0] for i in genes] # leave out numbers from first line
    header = "ST\t" + "\t".join(gene_names)
    #print(header)
    st_out.write(header+"\n")
    combined_out.write("Accession_nr\t"+header+"\n")
    mlst_infile.seek(0) #reset cursor at beginning of file

    new_st_count = 5000 #counting index for new sequence types

    # go through mlst_infile and extract info for theother two files
    for line in lines:
        accession, schema, st, gene1, gene2, gene3, gene4, gene5, gene6, gene7 = line.split('\t')
        
        # extract meta info
        accession = "GCA"+ accession.split("GCA")[1]
        accession = accession.split(".1.")[0]
        meta_out.write(accession+"\n")

        # extract genomic info
        genes = [gene1, gene2, gene3, gene4, gene5, gene6, gene7]
        allele_numbers = [i.split('(',1)[1] for i in genes] #remove front (
        allele_numbers = [i.split(')',1)[0] for i in allele_numbers] #remove back )

        if st == "-": #when the st is unknown by the mlst program
            all_numbers = str(new_st_count)+"\t"+"\t".join(allele_numbers)
            new_st_count +=1
        else:
            all_numbers = st+"\t"+"\t".join(allele_numbers)
        #print(all_numbers)

        #check for errors in the numbers
        if "," in all_numbers or "~" in all_numbers:
            pass
        else: 
            st_out.write(all_numbers+"\n")
            combined_out.write(accession+"\t"+all_numbers+"\n")
        #break
    print('output files generated: \n\tsequence_types.txt: contains all STs and their genotypes\n\tmetadata.txt: accession numbers of all STs in the same order as the sequence_types.txt file\n\tcombined_mlst.txt: both files combined with the accession numbers in the first column')
