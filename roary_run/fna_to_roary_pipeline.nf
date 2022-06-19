#!/usr/bin/env nextflow

fasta_files = Channel.fromPath(
    '/home/milena/megasync/masters_thesis/genome_sets_for_evaluation/closely_related_sections/cluster1/*.fna.gz')
/* extract fasta files from gbff*/

process FilterFasta{
    publishDir '/home/milena/megasync/masters_thesis/genome_sets_for_evaluation/closely_related_sections/cluster1'

    /*cpus 3*/

    input:
    file sequence from fasta_files

    output:
    file sequence into fasta_files_trimmed

    script:
    """
    filterContigs.pl 5000 $sequence overwrite
    """
}

/* do prokka annotation */
process ProkkaAnnotation{
    publishDir '/home/milena/megasync/masters_thesis/genome_sets_for_evaluation/closely_related_sections/cluster1'

    /*cpus 3*/

    input:
    file sequence from fasta_files_trimmed

    output:
    file "${sequence.baseName}.gff" into gff_files

    script:
    """
    gunzip -k -f $sequence
    prokka ${sequence.baseName}
    cd PROKKA_*
    mv *.gff ../"${sequence.baseName}.gff"
    cd ..
    rm -r PROKKA_*
    rm ${sequence.baseName}
    """
    /*Prokka saves all output files into a directory called PROKKA_[currentDate],
    I want to extract them into the publishDir*/
}

/* assemble a pan genome with roary */

process RoaryPangenome{
    publishDir '/home/milena/megasync/masters_thesis/genome_sets_for_evaluation/closely_related_sections/cluster1'
    
    errorStrategy 'ignore'

    input:
    file all_gff_files from gff_files.collect()

    output:
    file "*.Rtab"
    file "a*"
    file "b*"
    file "c*"
    file "summary_statistics.txt"
    file "*.fa*"
    file "pan_genome_results"

    script:
    """
    roary -p 16 ${all_gff_files}
    """
    /*
    roary: generate the pan genome

    potential second line here: 
    query_pan_genome -a union *.gff 
    query_pan_genome: make a file with all the clusters that roary has generated. 
        Each cluster is a line with either the gene name or the group name in the first part, 
        and then a tab separated list of all proteins in the cluster
    */
}

/*roary_pangenome.view()*/