import subprocess as sp
import os
import math
from xmlrpc.client import ProtocolError
import matplotlib.pyplot as plt
import matplotlib.colors as cols
import numpy as np
import statistics as stats
from tqdm import tqdm

# blast the protein sequences in question against the VFDB (virulene factor data base)

#makeblastdb 
# -in dataset.fasta 
# -dbtype nucl 
# -parse_seqids

#tblastn
# -db assembly.fasta
# -query protein_sequence.fasta
# > blast_assembly_sequence.txt

vfdb = "/home/milena/megasync/masters_thesis/VFDB_database/VFDB_setA_nt.fas" # nt = nucleotides
#vfdb = "/home/milena/megasync/masters_thesis/VFDB_database/VFDB_setB_pro.fas" # pro = protein sequences


infile_directory = "/home/milena/megasync/masters_thesis/genome_sets_for_evaluation/single_STs/ST_50_results/analyze_jumps"

################ make VFDB file into database ################

#sp.check_output('makeblastdb -in '+vfdb+ ' -dbtype nucl -parse_seqids', shell = True)


################ blast sequences against VFDB ################

os.chdir(infile_directory)
print(os.getcwd())
file_list = os.listdir()
contig_files = [i for i in file_list if '_contig_numbers_sorted.tsv' in i]
print(len(contig_files), 'contig files')

if False: # this takes a long time, only run once
    for file in contig_files:
        accession = file[0:13]
        print("==================== " + accession + " ====================")
        contig_file = open(file, "r")
        lines = contig_file.readlines()
        sequences_directory = accession+"_sequences"
        os.chdir(sequences_directory) # go into each directory that contains the protein sequences needed for blast
        for line in lines: # each contig no and their genes
            contig_no = line.split("\t")[0]
            print("contig number: "+contig_no)
            genes = line.split("\t")[1]
            gene_list = genes.split(",")
            for gene in tqdm(gene_list): # for each gene
                gene = gene.split("\n")[0]
                sequence_filename = accession+"_sequences_"+gene+".fa" # extract the sequence file
                blast_outfile = accession+"_"+gene+"_blast_against_VFDB.txt"
                sp.check_output('tblastn -db '+vfdb+ ' -query '+sequence_filename+' > '+blast_outfile, shell = True)
        os.chdir("..")

################ extract annotations from blast outfiles ################
# hits that have an e-value <0.1 will be counted as matches, but i will extract all

accessions = [i[0:13] for i in contig_files]


for i,accession in enumerate(accessions):
    print("==================== " + accession + " ====================", str(i)+"/"+str(len(accessions)))
    sequences_directory = accession+"_sequences"
    if os.path.isfile(accession+"_contig_numbers_sorted.tsv") == False:
        print(accession+"_contig_numbers_sorted.tsv skipped, does not exist")
        continue
    genes_in_contigs_file = open(accession+"_contig_numbers_sorted.tsv", "r")

    print("output file: "+accession+"_gene_jumps_annotations.tsv")
    outfile_annotations = open(accession+"_gene_jumps_annotations.tsv", "w")
    outfile_annotations.write("\"contig_no\"\t\"gene-family\"\t\"e-value\"\t\"general_category\"\t\"vfdb_annotation_string\"\n")
    
    # go through the genes in each contig
    lines = genes_in_contigs_file.readlines()
    for line in lines:
        contig_no = line.split("\t")[0]
        print("\tcontig: "+contig_no)
        genes = line.split("\t")[1]
        gene_list = genes.split(",")
        for gene in gene_list:
            vfdb_blast_filepath = sequences_directory+"/"+accession+"_"+gene+"_blast_against_VFDB.txt"
            if "\n" in vfdb_blast_filepath:
                vfdb_blast_filepath = vfdb_blast_filepath.replace("\n", "")
            #print(vfdb_blast_filepath)
            vfdb_blast_file = open(vfdb_blast_filepath, "r")
            vfdb_lines = vfdb_blast_file.readlines()
            # find indexes of string annotation lines (start with ">")
            # sometimes the annotation goes over multiple lines, so I will extract the two folloing ones as well
            annotation_exists = [line for line in vfdb_lines if ">" in line] #check if there are annotation lines
            if len(annotation_exists) == 0:
                annotation_lines = "no match found"
            else:
                start_annotation_line_index = [i for i, line in enumerate(vfdb_lines[0:100]) if ">" in line]
                end_annotation_line_index = [i for i, line in enumerate(vfdb_lines[0:100]) if "Length" in line] #length of the match is in the line after the annotation
                #print(start_annotation_line_index, end_annotation_line_index)
                if len(start_annotation_line_index) < 1:
                    annotation_lines = "no match found"
                else:
                    annotation_lines = vfdb_lines[start_annotation_line_index[0]:end_annotation_line_index[1]]
                    annotation_lines = [line.replace("\n", "") for line in annotation_lines]
                    annotation_lines = " ".join(annotation_lines).replace(">","")
                    #print(annotation_lines.split("-")) #the first [0] length would be the the length of the query
                    if len(annotation_lines.split("- ")) >1:
                        general_annotation = annotation_lines.split("- ")[1]
                        general_annotation = general_annotation.split(" (")[0]
                        general_annotation = general_annotation.strip()
                        general_annotation = general_annotation.replace("  ", " ")

            #line 22 contains the highest match -> index 21
            e_value = vfdb_lines[21].strip().split(" ")[-1] #strip() removes leading and trailing whitespace
            e_value = e_value.split("\n")[0]
            #print("E-value: "+e_value)

            line_in_outfile = "\""+contig_no+"\"\t\""+gene+"\"\t\""+e_value+"\"\t\""+general_annotation+"\"\t\""+annotation_lines+"\"\n"
            outfile_annotations.write(line_in_outfile)

            vfdb_blast_file.close()
    
    sequences_directory = accession+"_sequences"
    #os.chdir(sequences_directory) # go into each directory that contains the protein sequences needed for blast

    outfile_annotations.close()
    genes_in_contigs_file.close()