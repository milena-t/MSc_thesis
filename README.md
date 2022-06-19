# Large-scale pan-genome analysis on _Campylobacter jejuni_

The thesis that uses this code was written at Uppsala Universitet and can be found on [DiVA](http://uu.diva-portal.org/smash/record.jsf?dswid=-6916&pid=diva2%3A1671237&c=1&searchType=SIMPLE&language=en&query=large-scale+pan-genome+analysis+of+campylobacter+jejuni&af=%5B%5D&aq=%5B%5B%5D%5D&aq2=%5B%5B%5D%5D&aqe=%5B%5D&noOfRows=50&sortOrder=author_sort_asc&sortOrder2=title_sort_asc&onlyFullText=false&sf=undergraduate).
The dataset consists of all _Campylobacter jejuni_ genomes that were present in the GenBank database in Autumn 2021.


## Creating a phylogenetic tree

### Multi-locus sequence typing
The software that was used for MLST can be found on [GitHub](https://github.com/tseemann/mlst). It can take multiple input files either with * or as a space-separated list. The output was used as the default .tsv file, but with an .out file extension. The campylobacter schema was used (`--scheme campylobacter`). The program set a maximum input limit, so the input files were run in sections of 1000, and then merged using the MLST/merge_multiple_mlst_outputs.py script. This process also filters out genomes with uncertainty in some allele numbers. (The script MLST/prepare_mlst_for_phyloviz.py also does this for a single script without merging).
When there are unknown STs in the dataset, The script will also assign them a new ST number to these, starting at 50,000 to avoid overlap with existing ST numbers.

### Turn allele types into a newick tree 
the grapetree software ([website](https://achtman-lab.github.io/GrapeTree/), but it is a little out of date, see grapetree -h for up-to-date documentation) can be used to create a newick tree from a MLST output. It is important to mention that the headline of the file needs to start with '#', which is _not_ included in the merging script and has to be added manually! I used the `--profile` option.

## Creating a pan-genome

### Nextflow scripts to run roary

For the creation of the pan-genome, several processing steps were necessary before [roary](https://sanger-pathogens.github.io/Roary/) could be run, since roary requires a .gff file as input, which were annotated with [prokka](https://github.com/tseemann/prokka). The majority of the process, starting from fasta, was written in [nextflow](https://www.nextflow.io/). There is also a small bash script to extract fasta files from gbff. The general process is this:
1. (potentially) convert gbff to fasta with the [bioconvert](https://bioconvert.readthedocs.io/en/master/installation.html) mode [genbank2fasta](https://bioconvert.readthedocs.io/en/master/ref_converters.html#bioconvert.genbank2fasta.GENBANK2FASTA)
2. Run the nextflow script, which includes:
2.1 Filter contigs shorter than 5000 bp to get rid of potential contaminations (script `filterContigs.pl`)
2.2 Perform a prokka annotation on the resulting fasta files (keeping only the .gbff file for each sequence)
2.3 Run roary on all prokka output files

### Process the roary output





