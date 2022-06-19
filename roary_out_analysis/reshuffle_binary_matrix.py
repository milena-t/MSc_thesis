import sys
import numpy as np
import os
import subprocess as sp
import csv

binary_file = '/home/milena/megasync/masters_thesis/genome_sets_for_evaluation/single_STs/ST_48_results/bootstrap_matrices/gene_presence_absence.Rtab'


# number of bootstrap matrices to create
n = 10

for i in range(1,n):
    out_directory = binary_file.split('.Rtab')[0]+'bootstrap_'+str(i)+'.Rtab'

    binary_array = np.loadtxt(binary_file, delimiter = '\t', dtype=str)
    index_columns = np.arange(len(binary_array[0]))

    # random.shuffle manipulates the existing matrix, doesn't return a new variable
    np.random.shuffle(index_columns)
    index_columns = list(index_columns)

    # swap the first column index back to the beginning
    nr_at_one = index_columns[0]
    ind_of_first =  index_columns.index(0)
    index_columns[ind_of_first] = nr_at_one
    index_columns[0] = 0

    # reshuffle binary_array by column
    binary_array_remixed = binary_array[:,index_columns]

    # save array
    #binary_array_remixed.tofile(out_directory, sep = '\t')
    binary_list = list(binary_array_remixed)
    with open(out_directory, 'w+') as binary_out_file:
        for i, line in enumerate(binary_list):
            binary_out_file.write('\t'.join(line)+'\n')
    print('new binary matrix: '+out_directory)