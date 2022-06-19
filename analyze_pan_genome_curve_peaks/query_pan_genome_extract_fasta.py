import subprocess as sp
import os
import math
import matplotlib.pyplot as plt
import matplotlib.colors as cols
import numpy as np
import statistics as stats
import sys


#query_pan_genome 
# -g clustered_proteins
# -a gene_multifasta 
# -o genomes_that_cause_peaks
# -n group_3356,group_3102,group_7151,group_3870,group_3513,group_7727,group_1620,group_8231,group_11018,group_12204,group_5468,group_9539,group_12203,group_2744,group_10383,group_11017,group_564,rfbJ_2,group_10382,group_11011,group_11013,group_11015,hindIIIM,group_12186,group_12187,group_12198,group_12199,group_12200,group_12202,group_523,group_1419,group_3342,group_11016,group_1798,fmt_3,group_640,group_11019,group_3908,group_5481,group_8937,group_9578,group_10376,group_10377,group_10380,group_12191,group_12195,group_10853,group_12197,group_12201,group_12205,group_1611,group_2742,group_4290,group_4291,group_5555,group_6228,group_6248,group_8233,group_11012,group_11014,hindIIIM_2,group_12188,group_12189,group_12190,group_12192,group_12193,group_12194,group_12196,group_2215,group_3343,group_4355,group_522,group_9580 
# /home/milena/masters_thesis/code/pangenome_curve_peaks/genomes_that_cause_peaks/GCA_004913455.1.gbff.gff

gene_cluster_lists_directory = '/home/milena/megasync/masters_thesis/genome_sets_for_evaluation/single_STs/ST_48_results/analyze_jumps' 
# contains csv lists of relevant genes for each genome
clustered_proteins = '/home/milena/megasync/masters_thesis/genome_sets_for_evaluation/single_STs/ST_48_results/roary_output/clustered_proteins' 

genome_gff_directory = '/home/milena/megasync/masters_thesis/genome_sets_for_evaluation/single_STs/ST_48'

# python3 find_genes_in_peaks_in_the_pangenome_curve.py /home/milena/masters_thesis/code/pangenome_curve_peaks/genomes_that_add_too_many_genes_my_script.txt /home/milena/masters_thesis/jejuni_1000_-i90/gene_presence_absence.Rtab


os.chdir(gene_cluster_lists_directory)
file_list = os.listdir(gene_cluster_lists_directory)
gene_lists = [i for i in file_list if '.txt' in i and 'GCA' in i]
#print(gene_lists)
genomes_dir = os.listdir(genome_gff_directory)
genomes_list = [i for i in genomes_dir if '.gff' in i]
#print(genomes_list[0:20])

#loop through over gene csv
for gene_list in gene_lists:
    accession = gene_list[0:13]
    print(accession)
    #genome = [gene_cluster_lists_directory+'/'+i for i in genomes_list if accession in i][0]
    genome = [i for i in genomes_list if accession in i][0]
    #print(genome)
    #print(genome.split('/')[-1])
    sequences_out_folder = gene_list[0:13]+'_sequences' # add gene_cluster_lists_directory+'/'+ in front of the gene_list for absolute file path
    if os.path.isdir(sequences_out_folder) == False:
        os.mkdir(sequences_out_folder) # make output folder

    gene_list_file = open(gene_list, 'r') #read csv file with gene names
    gene_list = gene_list_file.readlines()[0]
    gene_list = gene_list[:-1] #remove trailing comma
    #print(gene_list)

    command_part1 = 'query_pan_genome -g '+clustered_proteins+' -a gene_multifasta'
    #command_part1 = 'query_pan_genome -a gene_multifasta'
    #print(command_part1)
    command_out_directory = ' -o '+sequences_out_folder# + '/' +accession
    #print(command_out_directory)
    #command_gene_csv = ' -n '+ genome_gff_directory+'/ '+gene_list + ' '
    command_gene_csv = ' -n '+gene_list + ' '
    command_gene_csv = command_gene_csv.replace("(", "")
    command_gene_csv = command_gene_csv.replace(")", "")
    #print(command_gene_csv)
    command_genome = genome_gff_directory+'/'+genome
    #print(command_genome)

    os.chdir(sequences_out_folder)
    print(os.getcwd())
    print(command_part1+command_out_directory+command_gene_csv+command_genome)
    sp.check_output(command_part1+command_out_directory+command_gene_csv+command_genome, shell=True)
    print('=============', accession, 'done =============')
    os.chdir('..')
   # break
