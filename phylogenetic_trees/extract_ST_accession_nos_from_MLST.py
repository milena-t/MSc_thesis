import subprocess as sp
import os
import threading
from types import new_class
import numpy as np
import sys

# take a list of STs that make up the clusters as input
# (each cluster one line, STs separated by commas. No comma at the end of a line!)
# input file as command line argument
# 
st_clusters_file = "/home/milena/megasync/masters_thesis/code/phylogenetic_trees/outbreaks_STs_3_biggest_clusters.txt"
st_clusters_file = "/home/milena/megasync/masters_thesis/code/phylogenetic_trees/outbreak_ST_48.txt"
st_clusters_file = "/home/milena/megasync/masters_thesis/code/phylogenetic_trees/long_sections_STs.txt"
# = sys.argv[1] 

# The second argument is the third output file from 'merge_mlst_outputs.py',
#  in which the sequence types and accesstion number are in one file (accessions in the first column)
# 
st_and_accessions_file = "/home/milena/megasync/masters_thesis/mlst/jejuni_all_analysis/combined_mlst.txt"
#  = sys.argv[2]

# get working directory (the one with st_clusters_file)
directory_levels = st_clusters_file.count('/')
working_directory = st_clusters_file.split('/')
working_directory = '/'.join(working_directory[:-1])
os.chdir(working_directory)
print('The working directory is:', os.getcwd())

with open(st_clusters_file, 'r') as st_clusters, open(st_and_accessions_file, 'r') as st_and_accessions:
	clusters = st_clusters.readlines()
	# remove newline character at the end of each cluster
	clusters = [i.split('\n')[0] for i in clusters]
	no_clusters = len(clusters)
	print('There are', no_clusters, 'clusters\n')
	
	# get lines from ST list
	st_and_accessions_list = st_and_accessions.readlines()[1:] # remove header
	# extract STs and ATs and acc. numbers
	st_list = [int(i.split('\t')[1]) for i in st_and_accessions_list]
	at_list = [i.split('\t', 2)[2] for i in st_and_accessions_list]
	acc_list = [i.split('\t')[0] for i in st_and_accessions_list]
	# merge all cluster accession numbers in a csv file for R
	R_csv = open('all_clusters_R.csv', 'w')
	
	#go through each cluster
	for i, cluster in enumerate(clusters): # i is the cluster index
		print('\n============== Cluster', i+1, 'is being processed ==============')
		# parse cluster to list of integers
		cluster = [int(i) for i in cluster.split(',')]
		print('There are', len(cluster), 'sequence types in the cluster')
		# find the indexes of each cluster element in the list of st types
		st_indices = []
		for cluster_st in cluster:
			#print([at_list[i] for i, st in enumerate(st_list) if st == cluster_st])
			st_indices.extend([acc_list[i] for i, st in enumerate(st_list) if st == cluster_st])
		with open('cluster'+str(i+1)+'.gclist', 'w') as out_accessions:
			print('There are', len(st_indices), 'Genomes in the cluster')
			all_accs_string_gclist = '\n'.join(st_indices)
			out_accessions.write(all_accs_string_gclist)	
		# append R csv file
		# The csv file is organized by rows, with the first element being the cluster label
		R_csv.write(','.join(st_indices)+'\n')
	R_csv.close()
	print("\t-> Outfiles here: "+working_directory+"/")