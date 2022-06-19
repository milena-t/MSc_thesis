import subprocess as sp
import os
import math
import matplotlib.pyplot as plt
import matplotlib.colors as cols
import numpy as np
import statistics as stats
import sys

from numpy.core.fromnumeric import mean

# /home/milenatrabert/masters_thesis/mapping_testset/fastq_10
# /media/milenatrabert/One\ Touch/rand1000fastqs/1_results/red_files
# directory as command line argument
# the arguments are organized like: python mapping.py[0] directory[1]

code_dir = os.getcwd()

directory = sys.argv[1] # [0] would be the name of this file
os.chdir(directory)
print(os.getcwd())

#extract all relevant files
all_files = os.listdir()
coverage_files = [i for i in all_files if '.coverage' in i and '.swp' not in i] # coverage files

print('nr. coverage files: ', str(len(coverage_files)), '\t', coverage_files[0], ' etc.') 

#summarize arrays of coverages and lengths into lists
contig_lengths = [] 
contig_coverages = []
#save names in correct order also into a list
contig_accessions = []

for file in coverage_files:
    accession, ending = file.split('.1', 1)
    print(file+' ============================')

    cov_file =  open(file, 'r')

    nr_contigs = len([0 for _ in cov_file]) -1 # exclude header line
    print('number of contigs:', nr_contigs)
    cov_file.seek(0) #reset cursor at beginning of file
    lines = cov_file.readlines()
    #lines = csv.reader(cov_file, delimiter='\t')
    #print(lines)
    index = 0
    contig_length = np.zeros(nr_contigs)
    contig_coverage = np.zeros(nr_contigs)
    lines = iter(lines) #convert to an iterator to use next()
    for line in lines:
        #print(line)
        if index == 0:
            next(lines) #skip header line
        else:
            line_elements = line.strip().split("\t")
            contig_length[index-1] = line_elements[2]
            contig_coverage[index-1] = line_elements[6]
        index += 1

    median_coverage = int(round(stats.median(contig_coverage)))
    print('empirical coverage:', median_coverage, '\n')
    contig_lengths.append(contig_length)
    contig_coverages.append(contig_coverage)
    contig_accessions.append(accession)
    cov_file.close()

### PLOT
color_list = []
if len(coverage_files)<11:
    color_list = [i for i in cols.TABLEAU_COLORS.values()] #list of 10 color values
else:
    #color_list = [i for i in cols.CSS4_COLORS.values()] #list of 148 color values
    color_list = ['#326496' for i in range(len(coverage_files))] 
    #color_list[-3:-1] = ['#A14A6B' for i in range(3)]

#set subplot dimensions
cols=4
rows=math.ceil(len(coverage_files)/cols)
plt.figure(figsize=(27, 12))
plt.subplots_adjust(wspace=1.0, hspace=0.5)
plt.suptitle("C: contig coverage vs. size, min 10 000 bp length", fontsize=18, y=0.95)


for i, accessions in enumerate(contig_accessions):
    ax = plt.subplot(rows,cols,i+1)
    plt.xlim(0, 400)
    ax.scatter(contig_coverages[i], contig_lengths[i], c = color_list[i])
    median_coverage = stats.median(contig_coverages[i])
    plt.axvline(x=median_coverage, c = color_list[i])

    ax.set_title(contig_accessions[i], c = color_list[i])
    #ax.text(185,max(contig_lengths[i]), 'coverage: '+str(int(round(median_coverage))), c = color_list[i])

#plt.supxlabel('coverage')
#plt.supylabel('contig length in bp')
plt.xlabel("coverage")
plt.ylabel("contig length in bp")
plt.show()


