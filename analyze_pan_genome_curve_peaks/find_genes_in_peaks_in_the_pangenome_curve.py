import subprocess as sp
import os
import math
import matplotlib.pyplot as plt
import matplotlib.colors as cols
import numpy as np
import statistics as stats
import sys


# take a list of files (file extension included) as input
# all in one line, comma separated
# /home/milena/megasync/masters_thesis/genome_sets_for_evaluation/single_STs/ST_50_results/analyze_jumps/jump_genome_names.txt
genome_name_list = "/home/milena/megasync/masters_thesis/genome_sets_for_evaluation/single_STs/ST_48_results/analyze_jumps/jump_genome_names.txt"

# Binary file returned by roary
# /home/milena/megasync/masters_thesis/genome_sets_for_evaluation/single_STs/ST_50_results/roary_output-i50/gene_presence_absence.Rtab
binary_matrix_file = "/home/milena/megasync/masters_thesis/genome_sets_for_evaluation/single_STs/ST_48_results/roary_output/gene_presence_absence.Rtab"

#/home/milena/megasync/masters_thesis/genome_sets_for_evaluation/single_STs/ST_50_results/analyze_jumps
out_directory = genome_name_list.split('/')[:-1]
out_directory = '/'.join(out_directory)
print('output_directory:\n'+out_directory)

# python3 find_genes_in_peaks_in_the_pangenome_curve.py /home/milena/masters_thesis/code/pangenome_curve_peaks/genomes_that_add_too_many_genes_my_script.txt /home/milena/masters_thesis/jejuni_1000_-i90/gene_presence_absence.Rtab

with open(genome_name_list, 'r') as genome_name_list, open(binary_matrix_file, 'r') as binary_matrix_file:
	genomes_that_cause_peaks = genome_name_list.readlines()[0].replace('\n', '').replace('.gz', '').split(',')
	#print(genomes_that_cause_peaks)
		#extract filenames and remove the newline characters. split into list

	binary_matrix_genes = binary_matrix_file.readlines()
	binary_matrix_genomes = binary_matrix_genes[0].replace('\n', '').split('\t')
		#extract genome names from binary matrix (1-indexed because the first element is the string 'Gene', which is the header for the gene names column)

	# find indices of all relevant genomes in the binary matrix
	genome_tics_dict = {genome : 1 for genome in genomes_that_cause_peaks}
		#dictionary with key(genome) : value(index)
	for genome in genome_tics_dict:
		genome_tics_dict[genome] = binary_matrix_genomes.index(genome)
	print('genomes that cause peaks: \n',genome_tics_dict)

	# go through each gene in these genomes:
	for genome in genome_tics_dict:
		genome_outfile = open(out_directory+'/'+genome[0:13]+'_genes_added_to_pangenome.txt', 'w')
		#genome_outfile.write('gene_name\tline_index\n')
		for line, gene in enumerate(binary_matrix_genes[1:]):
			gene = gene.replace('\n', '').split('\t')
			#check if the gene is in the genome
			if gene[genome_tics_dict[genome]] == '1':  
				if '1' not in gene[1:genome_tics_dict[genome]-1]:
					#print(gene[0], line)
					#genome_outfile.write(gene[0]+'\t'+str(line)+'\n')
					genome_outfile.write(gene[0]+',')
		print(genome[0:13]+'_genes_added_to_pangenome.txt', 'done')
		genome_outfile.close()
	print('all done!')
	