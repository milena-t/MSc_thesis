from pickle import FALSE
import subprocess as sp
import os
import math
from xmlrpc.client import ProtocolError
import matplotlib.pyplot as plt
import matplotlib.colors as cols
import numpy as np
import statistics as stats
from tqdm import tqdm
import re

# blast the protein sequences in question against the VFDB (virulene factor data base)

#makeblastdb 
# -in dataset.fasta 
# -dbtype nucl 
# -parse_seqids

#tblastn
# -db assembly.fasta
# -query protein_sequence.fasta
# > blast_assembly_sequence.txt

#target_dataset = "/home/milena/megasync/masters_thesis/genome_sets_for_evaluation/long_sections/cluster1/cluster1_multifasta.fna"
target_dataset = "/home/milena/megasync/masters_thesis/campy_plasmids_ncbi/campy_plasmids_multifasta.fna"
infile_directory = "/home/milena/megasync/masters_thesis/genome_sets_for_evaluation/long_sections/cluster1_results/analyze_jumps"


################ make target file into database ################

#sp.check_output('makeblastdb -in '+target_dataset+ ' -dbtype nucl -parse_seqids', shell = True)


################ blast sequences against database ################

os.chdir(infile_directory)
print(os.getcwd())
file_list = os.listdir()
contig_files = [i for i in file_list if '.fna' in i and "_contig_" in i and "plasmid_sequence" not in i]
print(len(contig_files), 'contig files')

if True: # this takes a long time, only run once
    for file in tqdm(contig_files):
        acc_contig_nr = file.split(".fna")[0]
        print("==================== " + file + " ====================")
        blast_outfile = acc_contig_nr+"_blast_against_ncbi_plasmids.txt"
        print("outfile: ", blast_outfile)
        sp.check_output('blastn -db '+target_dataset+ ' -query '+file+' > '+blast_outfile, shell = True)


################ extract plasmid hits from blast outfiles ################
# hits that have an e-value <e-4 will be counted as matches, but i will extract all

if True:

    file_list = os.listdir(infile_directory)
    acc_contig_files = [i for i in file_list if '_blast_against_ncbi_plasmids.txt' in i]
    print(len(contig_files), 'blast outfiles')

    outfile_plasmids = open("ncbi_plasmid_hits.tsv", "w")
    outfile_plasmids.write("acc_and_contig_no\tncbi_plasmid_name\tscore\te-value\n")

    for filename in acc_contig_files:

        print("==================== " + filename + " ====================")
        blast_file = open(filename, "r")
        blast_lines = blast_file.readlines()
        top_hit=""
        if len(blast_lines)<22:
            top_hit = "no_match                                           0\tx\n"
        else:
            top_hit = blast_lines[22]
            print(top_hit)
            if len(top_hit)<2:
                top_hit = "no_match                                           0\tx\n"
        outfile_plasmids.write(filename.split("_blast_against_ncbi_plasmids.txt")[0]+"\t"+top_hit)

    outfile_plasmids.close()