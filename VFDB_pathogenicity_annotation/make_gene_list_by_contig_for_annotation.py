import subprocess as sp
import os
import math
from xmlrpc.client import ProtocolError
import matplotlib.pyplot as plt
import matplotlib.colors as cols
import numpy as np
import statistics as stats
import sys

# takes the list returned by "blast_outfiles_extract_contig_numbers.py" and reformats it by 
# sorting by contig no (first column) and adding a csv list of gene names (second column)

infile_directory = "/home/milena/megasync/masters_thesis/genome_sets_for_evaluation/single_STs/ST_48_results/analyze_jumps"

os.chdir(infile_directory)
print(os.getcwd())
file_list = os.listdir()
contig_files = [i for i in file_list if '_contig_numbers_for_genes.tsv' in i]
print(len(contig_files), 'contig files ')


# go through each file
for contig_file in contig_files:
    print(contig_file)
    contig_file_old = open(contig_file, "r")
    accession = contig_file[0:13]
    print(accession)
    contig_sorted_file = open(accession+"_contig_numbers_sorted.tsv", "w")

    contig_dict = {}

    lines = contig_file_old.readlines()
    for line in lines[1:]:
        line = line.split("\n")[0]
        #print(line)
        gene = line.split("\t")[0]#.split("sequences_")[1] # needs to be added for the long sections, not necessary for single STs
        #print(gene)
        if len(line.split("\t")[1]) <2:
            continue
        else:
            #print(gene)
            contig = int(line.split("\t")[1])
        if contig not in contig_dict:
            contig_dict[contig] = [gene]
        contig_dict[contig].append(gene)
    #print(contig_dict)
    
    for contig_no in contig_dict:
        line = str(contig_no) + "\t" + ",".join(contig_dict[contig_no])
        contig_sorted_file.write(line+"\n")

    print(accession+"_contig_numbers_sorted.tsv", "done")
    contig_sorted_file.close()
    contig_file_old.close()
