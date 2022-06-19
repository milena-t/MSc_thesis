import subprocess as sp
import os
import math
import matplotlib.pyplot as plt
import matplotlib.colors as cols
import numpy as np
import statistics as stats
import sys

# Binary file returned by roary
# /media/milena/One\ Touch/jejuni_1000_-i90/gene_presence_absence.Rtab
binary_matrix_file = sys.argv[1]
out_suffix = sys.argv[2] #"cumulative_gene_content.tsv"
out_directory = binary_matrix_file.split('/')[0:-1]
out_directory = '/'.join(out_directory)
# python3 cumulative_size_of_pangenome_from_binary_matrix.py /media/milena/One\ Touch/jejuni_1000_-i90/gene_presence_absence.Rtab



with open(out_directory+'/'+out_suffix, 'w') as cumulative_gene_content_file, open(binary_matrix_file, 'r') as binary_matrix_file:
    cumulative_gene_content_file.write('genome\tpangenome_size\tgenomes_added\n')
    binary_matrix_file = binary_matrix_file.readlines()
    genomes = binary_matrix_file[0].replace('\n', '').split('\t')
    genomes_dict = {genome : 0 for genome in genomes[1:]}
    for line in binary_matrix_file[1:]:
        line = line.replace('\n', '').split('\t')
        first1 = line.index('1')
        genome_added = genomes[first1]
        genomes_dict[genome_added] +=1
    #print(genomes_dict)
    no_genomes_cumulative = np.cumsum(np.array(list(genomes_dict.values())))
    no_genomes_added = np.array(list(genomes_dict.values()))
    genome_names = list(genomes_dict.keys())
    #print(no_genomes_added)
    #print(no_genomes_cumulative)
    #print(genome_names)

    for i in range(len(genome_names)):
        cumulative_gene_content_file.write(genome_names[i]+'\t'+str(no_genomes_cumulative[i])+'\t'+str(no_genomes_added[i])+'\n')
print('output file: '+out_directory+'/'+out_suffix)