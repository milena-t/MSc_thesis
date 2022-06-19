import subprocess as sp
import os
import math
from xmlrpc.client import ProtocolError
import matplotlib.pyplot as plt
import matplotlib.colors as cols
import numpy as np
import statistics as stats
import sys


#makeblastdb 
# -in assembly.fasta 
# -dbtype nucl 
# -parse_seqids

#tblastn
# -db assembly.fasta
# -query protein_sequence.fasta
# > blast_assembly_sequence.txt

genomes_directory = '/home/milena/megasync/masters_thesis/genome_sets_for_evaluation/single_STs/ST_50_results/analyze_jumps' 

os.chdir(genomes_directory)
print(os.getcwd())
file_list = os.listdir()
sequence_directories = [i for i in file_list if '_sequences' in i]
print(len(sequence_directories), 'sequence directories')
genomes_list = [i for i in file_list if '.fna' in i and 'fna.n' not in i] # change fna <-> fasta if it doesn't work
print(len(genomes_list), 'genome fasta files')

for genome_file in genomes_list:
    accession = genome_file[0:13]
    # make a tsv table with gene names and in which contig they appear
    genes_and_contig_nos = open(accession+'_contig_numbers_for_genes.tsv', 'w')
    genes_and_contig_nos.write('protein_name\tcontig_number\n')

    print('================', genome_file, '================')
    
    sequence_directory = [i for i in os.listdir() if accession in i and 'sequences' in i][0]
    print('================', sequence_directory, '   ================')
    os.chdir(sequence_directory)

    blast_outfiles = os.listdir()
    blast_outfiles = [i for i in blast_outfiles if '_blast_outfile.txt' in i]
    
    # count files to calculate progress percentages later
    nr_files = len(blast_outfiles)
    nr_files_tenth = int(nr_files/10)
    print(' There are '+ str(nr_files)+ ' BLAST output files')

    # go through each blast outfile
    for i, blast_file in enumerate(blast_outfiles):
        protein_name = blast_file.split('_blast_outfile')[0]
        protein_name = protein_name.replace("sequences_", "")
        blast_file = open(blast_file, 'r')
        lines_blast_file = blast_file.readlines()
        if len(lines_blast_file)<20:
            genes_and_contig_nos.write(protein_name+'\t000\n')
        else:
            best_hit_contig = lines_blast_file[21] #21 is the index of the line with the first blast hit in the reference database
            best_hit_contig = best_hit_contig.split('1 ')[0][-4:-1] # maximum number of digits: 3 
            #print(best_hit_contig)

            genes_and_contig_nos.write(protein_name+'\t'+str(best_hit_contig)+'\n')

        blast_file.close()
        #break
        if i%nr_files_tenth==0:
            print(' ',round((i/nr_files)*100), '% completed')
    
    genes_and_contig_nos.close()
    #get back to the directory above the current sequence directory to start the loop over    
    os.chdir('..')
    #break
