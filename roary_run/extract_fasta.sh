#!/bin/sh

# command: 
# bash extract_fasta.sh /home/milena/megasync/masters_thesis/genome_sets_for_evaluation/long_sections/cluster2

gunzip $1/*.gz
for genome in $1/*.gbff
do 
    bioconvert genbank2fasta $genome "${genome%.gbff}.fna"
done
gzip $1/*.gbff
