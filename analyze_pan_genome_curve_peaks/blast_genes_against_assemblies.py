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

genomes_directory = '/home/milena/megasync/masters_thesis/genome_sets_for_evaluation/single_STs/ST_50' 
gene_cluster_lists_directory = '/home/milena/megasync/masters_thesis/genome_sets_for_evaluation/single_STs/ST_50_results/analyze_jumps'

os.chdir(gene_cluster_lists_directory)
print(os.getcwd())
file_list = os.listdir()
sequence_directories = [i for i in file_list if '_sequences' in i]
print(len(sequence_directories), 'sequence directories')
genomes_list = [i for i in os.listdir() if '.fna' in i and 'fna.gff' not in i]
print(len(genomes_list), 'genome fasta files')

#move fasta sequences to the analyze_jumps folder
for sequence_dir in sequence_directories:
    accession = sequence_dir[0:13]+'.1.fna.gz'
    sp.check_output('cp '+genomes_directory+"/"+accession+' .', shell = True)
    sp.check_output('gunzip -f '+accession, shell = True) # -f forces overwrite

print('--------- copied genomes -----------')

genomes_list = [i for i in os.listdir() if '.fna' in i and 'fna.gff' not in i]
print(len(genomes_list), 'genome fasta files')

for genome_file in genomes_list:
    accession = genome_file[0:13]

    print('================', genome_file, '================')
    
    
    # make blast database for genome
    sp.check_output('makeblastdb -in '+genome_file.replace(".gz", "")+ ' -dbtype nucl -parse_seqids', shell = True)

    sequence_directory = [i for i in os.listdir() if accession in i and 'sequences' in i][0]
    print('================', sequence_directory, '   ================')
    os.chdir(sequence_directory)
    protein_files = os.listdir()
    protein_files = [i for i in protein_files if '.fa' in i]
    nr_files = len(protein_files)
    nr_files_tenth = int(nr_files/10)
    print(' There are '+ str(nr_files)+ ' protein files')
    #protein_names = [i.split(accession)[2][1:-3] for i in protein_files]

    #make blast search for each protein
    for i, protein_file in enumerate(protein_files):
        protein_name = protein_file.split(accession)[-1][1:-3]
        #outfile = [i for i in protein_names if i in protein_file][0]+'_blast_outfile.txt'
        outfile = protein_name+'_blast_outfile.txt'
        #print('  '+outfile)
    
        command = 'tblastn -db ../' +genome_file+ ' -query ' +protein_file+ ' > ' +outfile
        sp.check_output(command, shell=True)
        #break
        if i%nr_files_tenth==0:
            print(' ',round((i/nr_files)*100), '% completed')
    #get back to the directory above the current sequence directory to start the loop over    
    os.chdir('..')
    #break
