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


infile_directory = "/home/milena/megasync/masters_thesis/genome_sets_for_evaluation/long_sections/cluster1_results/analyze_jumps"
os.chdir(infile_directory)
file_list = os.listdir()
blast_outfiles = [i for i in file_list if "_blast_against_source_dataset.txt" in i]
eval_file_name = "plasmid_search_eval.txt"
eval_file = open(eval_file_name, "w")
eval_file.write('plasmid_name\tno._of_occurences\taccessions_it_occurs_in\n')

for blast_outfile_name in blast_outfiles:
    blast_outfile = open(blast_outfile_name, "r")
    blast_outfile_lines = blast_outfile.readlines()
    eval_file_string = blast_outfile_name.split("_blast_against_source_dataset.txt")[0]+"\t"
    accessions = ""
    i = 0
    for i, line in enumerate(blast_outfile_lines[22:]):
        length = len(line)
        line = line.split("...")
        if length<2 : #break if the line is empty (end of list of matches is reached)
            break
        accession = [i[0:13] for i in line][0]#[:-1] #splitting by \t doesnt work? also :-1 to remove ' \n' at the end that results from this split
        accessions = accessions+accession+","
        values = [float(i) for i in line[1].split(" ") if len(i)>1]
        #print(i, accession, values)
        if values[1]>10e-5: #break if seq identity is too low
            break
    accessions = accessions[:-1]#remove last trailing comma
    eval_file_string = eval_file_string+str(i)+"\t"+accessions+"\n"
    eval_file.write(eval_file_string)
    print(eval_file_string[0:99])
    blast_outfile.close()
    #break