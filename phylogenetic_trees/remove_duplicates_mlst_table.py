from ete3 import Tree, TreeStyle
import os

# When plotting with ete3, a big problem is that the software explicitly plots each instance of a sequence type
# For this reason I will remove the duplicates in the ST file, which can then be transformed to a newick tree with grapetree

mlst_in_filepath = '/home/milena/mega/masters_thesis/mlst/jejuni_mlst_all_outfiles/mlst_all_ST.tsv'
outfile_filepath = '/home/milena/mega/masters_thesis/mlst/jejuni_mlst_all_outfiles/mlst_unique_ST.tsv'
st_duplications_filepath = '/home/milena/mega/masters_thesis/mlst/jejuni_mlst_all_outfiles/ST_duplications.tsv'

st_dict = {}

mlst_in = open(mlst_in_filepath, 'r')
mlst_unique_out = open(outfile_filepath, 'w')

tsv_file_lines = mlst_in.readlines()
mlst_unique_out.write(tsv_file_lines[0])
for line in tsv_file_lines[1:]: #skip header line
    st = line.split('\t')[0]# 0 for no accession number in front, 1 if there is an accession number
    if st in st_dict:
        st_dict[st] += 1
    else:
        st_dict[st] = 1
        mlst_unique_out.write(line)

mlst_unique_out.close()
mlst_in.close()

st_duplication_file = open(st_duplications_filepath, 'w')
st_duplication_file.write('ST\tnumber\n')
for key in st_dict:
    st_duplication_file.write(key+'\t'+str(st_dict[key])+'\n')
st_duplication_file.close()
