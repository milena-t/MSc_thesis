import subprocess as sp
import os
import threading
import sys

# /home/milenatrabert/masters_thesis/mapping_testset/fastq_10
# directory as command line argument
# the arguments are organized like: python3 mapping.py[0] directory[1]

code_dir = os.getcwd()

directory = sys.argv[1] # [0] would be the name of this file
os.chdir(directory)
print(os.getcwd())

#extract all relevant files
all_files = os.listdir()
raw_forward = [i for i in all_files if '_1.fastq.gz' in i] # raw reads of forward files
raw_reverse = [i for i in all_files if '_2.fastq.gz' in i] # raw reads of reverse files
reference = [i for i in all_files if 'fna.gz' in i]

print('nr. forward files: ', str(len(raw_forward)), '\t\t', raw_forward[0], ' etc.') 
print('nr. reverse files: ', str(len(raw_reverse)), '\t\t', raw_reverse[0], ' etc.')
print('nr. reference files: ', str(len(reference)), '\t', reference[0], ' etc.')

# define output vectors for each step
forward_alignment = []
reverse_alignment = []

# loop through each genome
for i in range(len(reference)):
	accession, ending = reference[i].split('.1', 1)
	print('\n============================ '+accession+' ============================')

	# the files in forward, reverse and accession are not in the same order, so I have to make sure
	# to use the correct raw read files for the reference genome
	forward_reads = [i for i in raw_forward if accession in i][0]
	print(forward_reads)
	reverse_reads = [i for i in raw_reverse if accession in i][0]
	print(reverse_reads)
	
	# indexing
	command = 'bwa index -p index '+reference[i]
	print(command+' ============================')
	sp.check_output(command, shell=True)
	
	# get sam file
	# out_sam = accession + '.1.sam'
	out_bam = accession + '.1.bam'
	command = 'bwa mem -t 16 -P index '+forward_reads+' '+reverse_reads+' | samtools view -u | samtools sort -o '+ out_bam
	#command = 'bwa sampe index '+out_forward+' '+out_reverse+' '+raw_forward[i]+' '+raw_reverse[i]+' | samtools view -b '+ out_bam
	sp.check_output(command, shell=True)
	print(out_bam+' ============================')

	# calculate coverage
	out_coverage = accession + '.1.coverage'
	command = 'samtools coverage '+out_bam + ' > ' + out_coverage
	sp.check_output(command, shell=True)
	print(out_coverage+' ============================')

os.chdir(code_dir)
