import subprocess as sp
import os
import math
from xmlrpc.client import ProtocolError
import matplotlib.pyplot as plt
import matplotlib.colors as cols
import numpy as np
import statistics as stats
from tqdm import tqdm

# merge all genes from ncbi into one multifasta file for the blast search 

infile_directory = "/home/milena/megasync/masters_thesis/genome_sets_for_evaluation/single_STs/ST_48"
outfile_name = "ST50_multifasta.fna"

os.chdir(infile_directory)
print(os.getcwd())
file_list = os.listdir()
plasmids = [i for i in file_list if '.1.fna' in i and "fna.gff" not in i]

outfile = open(outfile_name, "w") 
for plasmid_filename in plasmids:
    plasmid_name = plasmid_filename.split(".1.fna")[0]
    print(plasmid_name)
    with open(plasmid_filename, "r") as plasmid_file:
        plasmid_lines = plasmid_file.readlines()
        outfile.write(">"+plasmid_name+plasmid_lines[0].replace(">", " "))
        for line in plasmid_lines[1:]:
            if ">" in line:
                continue
            else:
                outfile.write(line)

outfile.close()