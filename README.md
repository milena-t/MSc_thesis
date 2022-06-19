# Large-scale pan-genome analysis on _Campylobacter jejuni_

The thesis that uses this code was written at Uppsala Universitet and can be found on [DiVA](http://uu.diva-portal.org/smash/record.jsf?dswid=-6916&pid=diva2%3A1671237&c=1&searchType=SIMPLE&language=en&query=large-scale+pan-genome+analysis+of+campylobacter+jejuni&af=%5B%5D&aq=%5B%5B%5D%5D&aq2=%5B%5B%5D%5D&aqe=%5B%5D&noOfRows=50&sortOrder=author_sort_asc&sortOrder2=title_sort_asc&onlyFullText=false&sf=undergraduate).
The dataset consists of all _Campylobacter jejuni_ genomes that were present in the GenBank database in Autumn 2021.


## Creating a phylogenetic tree

### Multi-locus sequence typing
The software that was used for MLST can be found on [GitHub](https://github.com/tseemann/mlst). It can take multiple input files either with * or as a space-separated list. The output was used as the default .tsv file, but with an .out file extension. The campylobacter schema was used (`--scheme campylobacter`). The program set a maximum input limit, so the input files were run in sections of 1000, and then merged using the MLST/merge_multiple_mlst_outputs.py script. This process also filters out genomes with uncertainty in some allele numbers. (The script MLST/prepare_mlst_for_phyloviz.py also does this for a single script without merging).
When there are unknown STs in the dataset, The script will also assign them a new ST number to these, starting at 50,000 to avoid overlap with existing ST numbers.

### Turn allele types into a newick tree 
the grapetree software ([website](https://achtman-lab.github.io/GrapeTree/), but it is a little out of date, see grapetree -h for up-to-date documentation) can be used to create a newick tree from a MLST output. It is important to mention that the headline of the file needs to start with '#', which is _not_ included in the merging script and has to be added manually! I used the `--profile` option.


## Handle the phylogenetic tree and extract datasets

All analyses were done with the [ETE3 toolkit](http://etetoolkit.org/). 

Many STs appear multiple times in the MLST output, which this is very inconvenient for plotting. Duplicates can be removed with the `phylogenetic_trees/remove_duplicates_mlst_table.py` script, and a new phylogenetic tree can be generated from only unique STs.
See my thesis for specifics on the datasets I decided to select. 

The selection results in lists of STs, which can be referred back to the MLST output to extract the accession numbers.


## Effects of filtering out small contigs
For a few select genomes, the raw reads were mapped back onto the contigs of the assembly to calculate a per-contig coverage, to potentially see contaminations, which are characterised by short contig length and low coverage. For more details see my thesis.


## Creating a pan-genome

### Nextflow scripts to run roary

For the creation of the pan-genome, several processing steps were necessary before [roary](https://sanger-pathogens.github.io/Roary/) could be run, since roary requires a .gff file as input, which were annotated with [prokka](https://github.com/tseemann/prokka). The majority of the process, starting from fasta, was written in [nextflow](https://www.nextflow.io/). There is also a small bash script to extract fasta files from gbff. The general process is this:
1. (potentially) convert gbff to fasta with the [bioconvert](https://bioconvert.readthedocs.io/en/master/installation.html) mode [genbank2fasta](https://bioconvert.readthedocs.io/en/master/ref_converters.html#bioconvert.genbank2fasta.GENBANK2FASTA)
2. Run the nextflow script, which includes:
  * Filter contigs shorter than 5000 bp to get rid of potential contaminations (script `filterContigs.pl`)
  * Perform a prokka annotation on the resulting fasta files (keeping only the .gbff file for each sequence)
  * Run roary on all prokka output files

### Process the roary output

All further analysis will be based on the roary output file "gene_presence_absence.Rtab", which is the binary matrix showing which gene families are present in which genomes. The pan-genome curve will be created in these steps:
1. potentially resuffle the binary matrix to create random permutations (pythons script)
2. get the cumulative size of of the pan-genome from the binary matrix (python script)
3. use milena-create_pangenome_plots.R script to plot the curve (modified from the existing roary script)

In this directory, there is also R files to create different plots included in my thesis for specific pan-genome outputs.


## Investigate pan-genome curve peaks

For details on this process, see my presentation, and especially slide 19 of the presentation. The process is this:
1. find which gene families are present in the jumps (`find_genes_in_peaks_in_the_pangenome_curve.py`)
2. Use the roary script "query_pan_genome" to extract the sequences of the gene families in the jumps (`query_pan_genome_extract_fasta.py`)
3. blast the gene sequences agains the respective genome in their jump (`blast_genes_against_assemblies.py`)
4. extract the contig number that each gene family matches (`blast_outfiles_extract_contig_numbers.py`)
5. process and plot the results (`contig_positions_of_added_genes_processing_and_plot.py`)


## Functional annotation

### VFDB

In the functional annotation directory, there are scripts to perform the pathogenicity annotation by blasting gene family sequences against a database of VFDB sequences. The previously extracted individual gene family sequences are the query that gets compared to the reference database of annotated VF sequences
1. `blast_against_VFDB_annotation.py` to sort the input files
2. `blast_against_VFDB_annotation` does what it says it does, also includes the extraction of the annotations from the blast out files
3. `plot_functional_annotation_contigs.R` plot the functional annotaion categories as a barplot


### NCBI plasmids

