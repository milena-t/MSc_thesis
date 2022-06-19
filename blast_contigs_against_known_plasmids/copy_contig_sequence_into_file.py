import subprocess as sp
import os
import math
from xmlrpc.client import ProtocolError
import matplotlib.pyplot as plt
import matplotlib.colors as cols
import numpy as np
import statistics as stats
from tqdm import tqdm

# use the contig_numbers_sorted.tsv file to get all the contigs that are relevant. 
# sum up all the genes, and only use contigs that have min 10% of all detected genes


infile_directory = "/home/milena/megasync/masters_thesis/genome_sets_for_evaluation/long_sections/cluster2_results/analyze_jumps"

os.chdir(infile_directory)
print(os.getcwd())
file_list = os.listdir()
contig_list = [i for i in file_list if 'contig_numbers_sorted.tsv' in i]

for contig_file_filename in contig_list:
    accession = contig_file_filename[0:13]
    print("\n================== "+ accession + " ==================")
    genome_filename = accession+".1.fasta" #change fna <-> fasta 
    with open(contig_file_filename, "r") as contig_file, open(genome_filename, "r") as genome_file:
        contig_nos = contig_file.readlines()
        contig_nos = [i.split("\t")[0] for i in contig_nos]
        fasta_lines = genome_file.readlines()
        # make dict with fasta headers and their line numbers/indices
        fasta_headers_dict = {i:line.split(".")[0]+" " for(i, line) in enumerate(fasta_lines)if ">" in line} #add " " so that -1 indexing works later
        #print(fasta_headers_dict)
        #print(contig_nos)
        for contig_no in contig_nos:
            output_fasta_name = accession+"_contig_"+contig_no+".fna"
            index_begin = 0
            index_end = 0
            for i in fasta_headers_dict:
                #print(fasta_headers_dict[i][-4:-1])
                if int(contig_no) == int(fasta_headers_dict[i][-4:-1]):
                    index_begin = i
                    #print(i)
                    print("from: ",fasta_headers_dict[i])
                if int(contig_no)+1 == int(fasta_headers_dict[i][-4:-1]):
                    index_end = i #not include header from following contig
                    #print(i)
                    print("to: ",fasta_headers_dict[i])
            if index_end == 0:
                index_end = len(fasta_lines)
            print("line indices in source file: ", index_begin, index_end)
            print("output filename:", output_fasta_name)      

            output_fasta = open(output_fasta_name, "w")
            #print("".join(fasta_lines[index_begin:index_begin+10]))
            output_fasta.write("".join(fasta_lines[index_begin:index_end]))
            output_fasta.close()
            print("\n")
    #break