import matplotlib.pyplot as plt
import matplotlib.colors as cols
from matplotlib import ticker
import numpy as np
import os
from math import *
import warnings
warnings.filterwarnings("ignore") # I AM EEEEEEVIL

def extract_contig_numbers(table_file_name):#, contig_lengths):
  #extract a dictionary of all the contig numbers of genes in 'table' and how many times they occur
  table_file = open(table_file_name, 'r')
  table_lines = table_file.readlines()
  if len(table_lines) == 0:
    return([])
  contigs = []
  for line in table_lines:
    contigs.append(line.split('\t')[1].split('\n')[0])
    #if int(line.split('\t')[1].split('\n')[0]) not in contig_lengths:
    #  print(line)
  contigs = contigs[1:]
  #print(len(contigs))
  contigs = np.array([int(i) for i in contigs if len(i)>0])
  #print(len(contigs))
  #print(contigs)
  contig_occurance = dict.fromkeys(contigs)
  contig_occurance = {x:0 for x in contig_occurance}
  for contig in contigs:
    contig_occurance[contig] += 1
  # keys: contig numbers
  # values: how often they appear
  table_file.close()
  return(contig_occurance)

def get_contig_size(fasta_file_name):
  fasta_file = open(fasta_file_name, 'r')
  fasta_lines = fasta_file.readlines()
  contig_lengths = {}
  line_no = 0
  contig_no = 0
  for line in fasta_lines:
    if '>' in line[0:5]: #only in the beginning of the line, not in the annotation string
      contig_no = int(line.split(' ')[0][-5:-2])
      contig_lengths[contig_no] = 0
      contig_lengths[contig_no-1] = line_no*60 # there are 60 bases in a line
      #print(contig_no-1, line_no)
      line_no = 0
    line_no +=1
  contig_lengths[contig_no] = line_no*60
  fasta_file.close()
  contig_lengths.pop(0)
  return(contig_lengths)

# get all contig tables and fasta files
code_dir = os.getcwd()
analyze_directory = '/home/milena/megasync/masters_thesis/genome_sets_for_evaluation/single_STs/ST_50_results'
#analyze_directory = '/home/milena/megasync/masters_thesis/genome_sets_for_evaluation/long_sections/cluster1_results'
tables_directory = analyze_directory+'/analyze_jumps'
os.chdir(tables_directory)
all_files = os.listdir(tables_directory)
contig_tables = [i for i  in all_files if '_contig_numbers_for_genes.tsv' in i]
fasta_files = [i for i in all_files if '.1.fna' in i and 'fna.n' not in i and 'fna.gff' not in i] # change fna <-> fasta if it doesn't work
print(fasta_files)
accessions = [i.split('_contig_numbers_for_genes.tsv')[0] for i in contig_tables]




if False:
  contig_genes = extract_contig_numbers('GCA_004995305_contig_numbers_for_genes.tsv')
  print(contig_genes)
  contig_lengths = get_contig_size(fasta_files[0])
  print(contig_lengths)

  #calculate how many genes are added per 10 000 bp(arbitrary choice, just to get some prettier values) for each contig
  for contig in contig_genes: 
    if contig not in contig_lengths:
      continue #skip if contig key not in contig lengths, something weird with my string parsing?
    contig_genes[contig]=contig_genes[contig]/contig_lengths[contig]*10000
  #print(contig_genes)



############## PLOT ##############

# if only certain accessions should be shown
if False:
  plots_to_show = [1,2,3,5,7,8,11,12]
  accessions = [accessions[i-1] for i in plots_to_show]

if True:

  #set subplot dimensions
  cols=6
  rows = ceil(len(accessions)/cols)
  fig_width = 3*cols
  fig_height = 3*rows
  #fig = plt.figure(figsize=(12, 6)) #for two rows
  #fig = plt.figure(figsize=(12, 12)) #for three rows
  fig = plt.figure(figsize=(fig_width, fig_height))
  plt.subplots_adjust(wspace=0.8, hspace=0.8)
  plt.suptitle("genes added to pangenome per contig", fontsize=20, y=0.92)

  # make big subplot for common axis labels
  common = fig.add_subplot(111)
  common.set_xlabel('Genes added per 10000 basepairs')
  common.set_ylabel('Contig Numbers')

  # make specific colors for the mis-annotated genomes
  # wrongly annotated genomes are:
  other_species = ['GCA_005078665', 'GCA_004995305']
  header_colors = ['#A0A0A0' if i in other_species else '#000000' for i in accessions]
  print(header_colors, len(header_colors))

  for i, accession in enumerate(accessions):
    print(accession)
    ax = plt.subplot(rows,cols,i+1)
    # get dictionaries for each accession
    contig_gene_file = [i for i in contig_tables if accession in i][0]
    #print(contig_gene_file)
    # all contig numbers for each genome
    contig_genes = extract_contig_numbers(contig_gene_file)
    if len(contig_genes) == 0:
      pass
    contig_length_file = [i for i in fasta_files if accession in i][0]
    print(contig_length_file)
    #print()
    contig_lengths = get_contig_size(contig_length_file)
    no_contigs = len(contig_lengths)
    for contig in contig_genes:
      if contig not in contig_lengths:
        continue #skip if contig key not in contig lengths, something weird with my string parsing?
      if contig_lengths[contig] == 0:
        continue
      contig_genes[contig]=contig_genes[contig]/contig_lengths[contig]*10000

    ax.bar([str(i) for i in contig_genes.keys()], contig_genes.values())#, c = color_list[i])

    ax.set_title(accession+', '+str(no_contigs)+' contigs', c = header_colors[i], fontsize = 12)

    #change x axis tics font size
    axis_ticks_fontsize = 11
      # x ticks
      # adjust smaller if there is a lot of ticks
    if len(contig_genes.keys())>8:
      axis_ticks_fontsize = 10
    if len(contig_genes.keys())>12:
      axis_ticks_fontsize = 9
    if len(contig_genes.keys())>14:
      axis_ticks_fontsize = 8
    if len(contig_genes.keys())>16:
      axis_ticks_fontsize = 6
    if len(contig_genes.keys())>18:
      axis_ticks_fontsize = 5
    ax.set_xticklabels([str(i) for i in contig_genes.keys()], fontsize=axis_ticks_fontsize)
      # y ticks 
    axis_ticks_fontsize = 14
    print(contig_genes)
    no_ticks = 4 # set number of ticks
    ticks = range(0, ceil(max(contig_genes.values()))*no_ticks, ceil(max(contig_genes.values())))
    ticks = [round(i/no_ticks, 1) for i in ticks]
    print(ticks)
    yticks = ticker.MaxNLocator(no_ticks) #set max number of axis ticks
    ax.yaxis.set_major_locator(yticks)
    ax.set_yticklabels(ticks, fontsize=axis_ticks_fontsize)
    if ticks[-1] < 1:
      ax.set_ylim([0,1.1])

  #plt.supxlabel('coverage')
  #plt.supylabel('contig length in bp')
  #plt.ylabel("Genes added per 10000 basepairs")
  #plt.xlabel("Contig Numbers")

  fig.text(0.5, 0.01, 'Contig Numbers', ha='center', fontsize = 14)
  fig.text(0.01, 0.5, 'Number of genes added per 10000 basepairs', va='center', rotation='vertical', fontsize = 14)

  os.chdir(code_dir)
  fig.tight_layout(h_pad=4, w_pad = 2) #adjust w_pad to -3 if a y label is added
  plt.subplots_adjust(top=0.85, left = 0.075, bottom = 0.1)
  fig.savefig(analyze_directory+'/genes_added_per_contig.png', format = 'png', dpi = 300)
  #plt.show()

  