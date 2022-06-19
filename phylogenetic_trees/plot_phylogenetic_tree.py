from ete3 import *
#import matplotlib.colors as mcolors

#import tree
tree_file_path = '/home/milena/megasync/masters_thesis/mlst/jejuni_all_analysis/mlst_unique_st_tree.new'
tree_file = open(tree_file_path, 'r')
tree_string = tree_file.readlines()[0].split('\n')[0]
tree = Tree(tree_string)

#import list of cluster STs
st_clusters = []
size_sd_clusters = []
with open("outbreaks_STs_smaller_clusters.txt", 'r') as outbreak_sts_file:
    outbreak_clusters = outbreak_sts_file.readlines()
    for cluster in outbreak_clusters:
        st_clusters.append(cluster[:-1].split(',')) #(:-1 to remove the \n char at the end of each line)
        size_sd_clusters.append(len(st_clusters[-1]))
with open("outbreaks_STs_3_biggest_clusters.txt", 'r') as outbreak_sts_file:
    outbreak_clusters = outbreak_sts_file.readlines()
    for cluster in outbreak_clusters:
        st_clusters.append(cluster[:-1].split(',')) #(:-1 to remove the \n char at the end of each line)
        size_sd_clusters.append(len(st_clusters[-1]))

def st_layout(node): #mark ST48 and ST50
    STs = set(["48", "50"])
    if node.name in STs:
        print(node)
        node.img_style["size"] = 200000
        node.img_style["fgcolor"] = "red"
    else:
        node.img_style["size"] = 0
        node.img_style["fgcolor"] = "white"


# plot the phylogenetic tree
if True:

    #cols = ['rosybrown', 'indianred', 'darksalmon', 'peru', 'navajowhite', 'goldenrod', 'darkseagreen', 'lightgreen', 'mediumaquamarine', 'skyblue', 'cornflowerblue', 'mediumpurple', 'thistle', 'palevioletred']
    cols_outbreak_clusters = ['orchid', 'thistle', 'palevioletred']
    cols_long_sections = ['lightsteelblue', 'cornflowerblue']
    col_st48 = ['darkseagreen']
    cols = cols_outbreak_clusters + cols_long_sections + col_st48
    c = 0
    #print(cols)
    #print(cols[c])

    #remove blue dots from nodes
    nstyle = NodeStyle()
    #nstyle["size"] = 0
    nstyle["vt_line_width"] = 200 #use 10 when trying to see the leaf names, otherwise at least 200
    nstyle["hz_line_width"] = 200
    nstyle["size"] = 2000
    nstyle["fgcolor"] = "white"
    for n in tree.traverse():
        if n.name in set(["48", "50"]):
            print(n)
            nstyle["size"] = 200000
            nstyle["fgcolor"] = "red"
            nstyle["node_bgcolor"] = "red"
            nstyle["partition_bgcolor"] = "red"
            nstyle["faces_bgcolor"] = "red"
            #n.set_style(nstyle)
        else:
            nstyle["size"] = 0
            nstyle["fgcolor"] = "white"
        n.set_style(nstyle)

    if False: #TODO doesn't work yet
        #st_48 = ['48']
        single_st_node = tree.search_nodes(name="50")[0] # 50 or 48
        #n = tree.get_common_ancestor(["67818", "50"])
        print(single_st_node)
        print(n)
        nst = NodeStyle()
        #nst3['bgcolor'] = cols[c]
        nst['fgcolor'] = cols[c]
        nst["size"] = 20000
        #nst3['size'] = 0
        nst["shape"] = "circle"
        single_st_node.set_style(nst)

    if False:
        single_st_node = tree.search_nodes(name="50")[0] # 50 or 48
        single_st_node.add_face(CircleFace(2000, "red", ), column = 1, position = "branch-right")   

    if False: #add markings for clusters of outbreaks
        for i in range(len(st_clusters)):
            if size_sd_clusters[i]>50:
                nst = NodeStyle()
                nst['bgcolor'] = cols[c]
                nst["size"] = 0
                n = tree.get_common_ancestor(st_clusters[i])
                n.set_style(nst)
                c = c+1

    if False: #add markings for clusters of long sections
        leaves_bottom_part = ['4349', '10086', '5271', '8373', '8361']
        nst1 = NodeStyle()
        nst1['bgcolor'] = cols[c]
        nst1["size"] = 0
        n = tree.get_common_ancestor(leaves_bottom_part)
        n.set_style(nst1)

        c = c+1

        #leaves_top_part = ['7843', '5790', '841']
        leaves_top_part = ['67818', '5673']
        nst2 = NodeStyle()
        nst2['bgcolor'] = cols[c]
        nst2["size"] = 0
        n = tree.get_common_ancestor(leaves_top_part)
        n.set_style(nst2)
    
        c = c+1 

    ts = TreeStyle()
    ts.show_leaf_name = False
    ts.mode = "c"
    ts.rotation = 90
    ts.arc_start = -180 # 0 degrees = 3 o'clock
    ts.arc_span = 350
    #ts.title.add_face(TextFace("Campylobacter jejuni", fsize=20), column=0)

    #tree.show(tree_style=ts)

    plot_name = "jejuni_MLST_tree_all_sections.png"
    tree.render(plot_name, h=10000, tree_style=ts) #width and height in pixels
    print('plot saved as',plot_name)
    