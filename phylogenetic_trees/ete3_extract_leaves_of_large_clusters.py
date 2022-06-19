from ete3 import Tree, TreeStyle, NodeStyle
import os

tree_file_path = '/home/milena/megasync/masters_thesis/mlst/jejuni_all_analysis/mlst_unique_st_tree.new'
#tree_file_path = '/home/milena/masters_thesis/mlst/jejuni1000_analysis/jejuni1000_grapetree_tree.new'
#tree_file_path = '/media/milenatrabert/One Touch/mlst/jejuni_mlst_all_outfiles/mlst_grapetree_tree.new'

tree_file = open(tree_file_path, 'r')
tree_string = tree_file.readlines()[0].split('\n')[0]
tree = Tree(tree_string)

#test tree
t = Tree("((((a1,a2,a3,a4)a5,a6)aa, (b1,b2)b3)ab, (c1, (d1,d2)d3)cd);", format=1)


# return clusters
#nodes, leaves, size = search_by_size(tree, size=i)

leaves_top_part = ['67818', '5673']
n2 = tree.get_common_ancestor(leaves_top_part)
all_leaves_top_part = n2.get_leaf_names()
print(len(all_leaves_top_part), 'leaves in the top part')
print(all_leaves_top_part, '\n')

leaves_bottom_part = ['4349', '10086', '5271', '8373', '8361']
n1 = tree.get_common_ancestor(leaves_bottom_part)
all_leaves_bottom_part = n1.get_leaf_names()
print(len(all_leaves_bottom_part), 'leaves in the bottom part:')
print(all_leaves_bottom_part) 


if True:
    # print leaf names to list that can be handled by the "extract_ST_accession_nos.py" script in the MLST folder
    # (each cluster one line, STs separated by commas. No comma at the end of a line!)
    with open("long_sections_STs.txt", 'w') as section_sts_file:
        for st_cluster in [all_leaves_top_part, all_leaves_bottom_part]:
            st_string = ','.join(st_cluster)
            #print(st_string)
            section_sts_file.write(st_string+'\n')
            #break

