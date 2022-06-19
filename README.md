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



