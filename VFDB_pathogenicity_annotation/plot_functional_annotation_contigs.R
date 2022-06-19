# make charts of funcitonal categories in contigs

#install.packages("viridis")     
library("viridis")
library("RPMG")
library("stringr")

################### FUNCTIONS FOR PROCESSING ###################

# extract all files in filepath that contain a substring (like ls *substring*)
get_file_with_substring <- function(substring, filepath){
  outfiles <- c()
  all_files <- list.files(filepath)
  for (file in all_files){
    if (grepl(substring, file, fixed = T)){
      outfiles[length(outfiles)+1] <- file
    }
  }
  return(outfiles)
}

# turn the *gene_jumps*annotations tables into a form that can be used with the baseR barplot
# the result is a table with the categories in the rows, and the contigs in the columns, and how many genes in each contig 
#   fit into each category.
# it contains several filtering steps
# (double checked with the ggplot function that can take the relevant columns of the raw table)
# order_vector makes sure that the rows of each output matrix are in the same order
make_stacked_barplot_table <- function(filename, order_vector){
  table <- read.table(filename, header = T)
  # process table
  table$general_category[table$e.value>10e-4] <- "no match found" #exclude matches with too high e-values
  table$contig_no <- as.factor(table$contig_no)
  table$general_category <- as.factor(table$general_category)
  # filter out contigs that have less than 10% of the total amount of genes (nrow(table)) in them
  plasmid_contigs <- strtoi(levels(table$contig_no)[summary(table$contig_no)>0.1*nrow(table)])
  if (length(plasmid_contigs)<1){ #when no single contig contains more than 10% of all genomes
    return("N")
  }
  
  # make barplot matrix to use baseR barplot function
  barplot_table <- matrix(0, length(levels(table$general_category)), length(plasmid_contigs))
  rownames(barplot_table) <- levels(table$general_category)
  colnames(barplot_table) <- plasmid_contigs
  # add amounts of annotations for each annot in each contig
  for (i in 1:length(plasmid_contigs)){
    values = summary(table$general_category[table$contig_no == plasmid_contigs[i]])
    barplot_table[,i] <- values
    # barplot_table[,i] <- values/sum(values) #(for relative values)
  }
  
  # rearrange matrix according to order_vector
  # make empty matrix with all annotations
  barplot_table_empty <- matrix(0, length(order_vector), length(plasmid_contigs))
  rownames(barplot_table_empty) <- order_vector
  colnames(barplot_table_empty) <- plasmid_contigs
  for (fun_annot in order_vector){
    if (any(grepl(fun_annot, rownames(barplot_table), fixed = T))){
      barplot_table_empty[fun_annot,] <- barplot_table[fun_annot,]
    }
  }
  barplot_table <- barplot_table_empty
  
  if (F){
    #  if the plot has all annotations in the order_vector
    if (nrow(barplot_table) == length(order_vector)){
      barplot_table <- barplot_table[order_vector,]
    }else{ # if some annotations don't appear
      # make empty matrix with all annotations
      barplot_table_empty <- matrix(0, length(order_vector), length(plasmid_contigs))
      rownames(barplot_table_empty) <- order_vector
      colnames(barplot_table_empty) <- plasmid_contigs
      for (fun_annot in order_vector){
        if (any(grepl(fun_annot, rownames(barplot_table), fixed = T))){
          barplot_table_empty[fun_annot,] <- barplot_table[fun_annot,]
        }
      }
      barplot_table <- barplot_table_empty
    }
  }
  
  
  return(barplot_table)
}


################### READING DATA ###################
"ncbi_plasmid_hits_only_hits.tsv"
directory <- "/home/milena/megasync/masters_thesis/genome_sets_for_evaluation/single_STs/ST_50_results/analyze_jumps"
directory <- "/home/milena/megasync/masters_thesis/genome_sets_for_evaluation/long_sections/cluster1_results/analyze_jumps"

order_vector <- c("Adherence",
                  "Biofilm",
                  "invasion",
                  "Effector delivery system",
                  "Exotoxin",
                  "Immune modulation",
                  "Motility",
                  "Nutritional/Metabolic factor",
                  "Stress survival",
                  "Exoenzyme",
                  "no match found")

extract_all_tables_list <- function(directory){

  setwd(directory)
  substring = "_gene_jumps_annotations.tsv"
  file_list <- get_file_with_substring(substring, directory)
  
  # extract a list of all barplot tables
  order_vector <- c("Adherence",
                    "Biofilm",
                    "invasion",
                    "Effector delivery system",
                    "Exotoxin",
                    "Immune modulation",
                    "Motility",
                    "Nutritional/Metabolic factor",
                    "Stress survival",
                    "Exoenzyme",
                    "no match found")
  
  all_tables_list <- list(0)
  
  #filename <- file_list[i] # for debugging purposes
  
  i=1
  for (file in file_list){
    # all_tables_list[[length(all_tables_list)+1]] <- make_stacked_barplot_table(file)
    print(file)
    print(i)
    barplot_table <- make_stacked_barplot_table(file, order_vector)
    print(barplot_table)
    if (length(barplot_table) > 1){
      all_tables_list[[i]] <- barplot_table
      names(all_tables_list[[i]]) <- substring(file, 1, 14)
    }
    i = i+1
  }
  
  len <- length(all_tables_list)
  return(all_tables_list)

}


############### SUM UP ALL CATEGORIES #############################

sum_up_categories <- function(all_tables_list){
  sums <- rep(0, length(order_vector))
  names(sums)<- order_vector
  for (table in all_tables_list){
    print(table)
    for(annot_category in order_vector){
      sums[annot_category] <- sums[annot_category]+sum(table[annot_category,])
    }
  }
  return(sums)
}

directory_st50 <- "/home/milena/megasync/masters_thesis/genome_sets_for_evaluation/single_STs/ST_50_results/analyze_jumps"
directory_st48 <- "/home/milena/megasync/masters_thesis/genome_sets_for_evaluation/single_STs/ST_48_results/analyze_jumps"
directory_g1 <- "/home/milena/megasync/masters_thesis/genome_sets_for_evaluation/long_sections/cluster1_results/analyze_jumps"
directory_g2 <- "/home/milena/megasync/masters_thesis/genome_sets_for_evaluation/long_sections/cluster2_results/analyze_jumps"

all_tables_list_st50 <- extract_all_tables_list(directory_st50)
filter_table <- read.table(paste(directory_st50, "/ncbi_plasmid_hits_only_hits.tsv", sep = "")) #keep only plasmids with hits
filter_table <- filter_table[[1]]

for (table in all_tables_list_st50){
  accession_table <- names(table)[[1]]
  print(paste("table", accession_table))
  contigs <- colnames(table)
  for (cont_nr in contigs){
    acc_cont <- paste(accession_table, cont_nr, sep = "contig_")
    print(acc_cont)
    index <- grep(acc_cont, filter_table, ignore.case = T)
    print(index)
  }
  if (length(index)<1){
    
  }
  
#  for(acc_present in filter_table){
#    accession_filter <- substr(acc_present, 1, 14)
#    #print(accession_filter)
#    if (accession_table == accession_filter){ #if the accession is represented
#      print(acc_present)
#
#    }
#  }
  #break
}

sums_st50 <- sum_up_categories(all_tables_list_st50)
norm_sum_st50 <- sums_st50/sum(sums_st50)

all_tables_list_st48 <- extract_all_tables_list(directory_st48)
sums_st48 <- sum_up_categories(all_tables_list_st48)
norm_sum_st48 <- sums_st48/sum(sums_st48)

all_tables_list_g1 <- extract_all_tables_list(directory_g1)
sums_g1 <- sum_up_categories(all_tables_list_g1)
norm_sum_g1 <- sums_g1/sum(sums_g1)

all_tables_list_g2 <- extract_all_tables_list(directory_g2)
sums_g2 <- sum_up_categories(all_tables_list_g2)
norm_sum_g2 <- sums_g2/sum(sums_g2)



###### plot ######
par(mfrow = c(2,3))
colors <- c(rainbow.colors(length(order_vector)-1), "gray")
# pie chart radius 
r <- 0.9
# main headings font size
msize = 1.5

barplot_sums <- function(norm_sum_st50,name){
  barplot(norm_sum_st50[-length(norm_sum_st50)], 
              col = colors[-length(norm_sum_st50)], 
              main = name,
              #sub = paste(as.character(round(sum(norm_sum_st50[-length(norm_sum_st50)])*100), 3), "% of all genes have a pathogenic annotatoin"),
              xlab = " ",  xaxt="n",
              ylab = "proportion",
              ylim = c(0, 0.085),
              cex.lab = 1.25, 
              cex.main = 1.5, cex.sub = 1.3)
        subtitle <- paste(as.character(round(sum(norm_sum_st50[-length(norm_sum_st50)])*100, 1)), "% of all genes have a pathogenic annotation")
        mtext(side=1, line=1, at=-0.3, adj=0, cex=0.9, subtitle)        
}

#pie(sums_g1, rep("", length(sums_g1)), radius = r, main = "Gene functions Group 1", cex.main = msize, col = c(rainbow.colors(length(order_vector)-1), "gray"))
barplot_sums(norm_sum_st50, "Proportion of functional \ncategories ST 50")
#pie(sums_g2, rep("", length(sums_g2)), radius = r, main = "Gene functions Group 2", cex.main = msize, col = c(rainbow.colors(length(order_vector)-1), "gray"))
barplot_sums(norm_sum_st48, "Proportion of functional \ncategories ST 48")

##### legend #####
#use legend from below plot sectoin

par(mar = c(0, 0, 1, 0)) #bottom,left,top,right
plot(1,1, type = "n", bty="n",
     mar = c(1,1,1,1),
     main = "Legend", xpd = TRUE, cex.main = msize,
     xaxt = "n", yaxt = "n",
     xlab = c(" "), ylab = c(" "))

legend(0.65,1.4,  inset=c(-0.1, 0),
       legend = order_vector[-length(norm_sum_st50)], 
       fill = colors[-length(norm_sum_st50)],
       cex = 1.5)
par(mar= c(5, 4, 4, 2) + 0.1) #reset margins to default
###### continue plot

#pie(sums_st48, rep("", length(sums_st48)), radius = r, main = "Gene functions ST 48", cex.main = msize, col = c(rainbow.colors(length(order_vector)-1), "gray"))
barplot_sums(norm_sum_g1, "Proportion of functional \ncategories group 1")
#pie(sums_st50, rep("", length(sums_st50)), radius = r, main = "Gene functions ST 50", cex.main = msize, col = c(rainbow.colors(length(order_vector)-1), "gray"))
barplot_sums(norm_sum_g2, "Proportion of functional \ncategories group 2")

############### PLOT INDIVIDUAL JUMPS ###############

# clear current plot first!
#dev.off(dev.list()["RStudioGD"])

#directory <- "/home/milena/megasync/masters_thesis/genome_sets_for_evaluation/single_STs/ST_48_results"
#png("functional_annotation_ST48.png")


colors = c(inferno(length(order_vector)-1), "gray")
colors = c(rainbow.colors(length(order_vector)-1), "gray")
len = as.integer(length(all_tables_list))+1 
#len = 23 #set for custom split (plot only 15 subplots)
# one plot
# par(mar=c(5.5, 4.5, 4, 14), xpd=TRUE) #make space for the legend to the right of the plot
# grid of multiple plots
cols = 4
rows = as.integer((len/cols)) +1 # maybe
#par(mfrow = c(cols,rows)) 
par(mfrow = c(rows,cols)) 
par(mar = c(5.1, 4.1, 4.1, 2.1))
#par(mar = c(1.1, 4.1, 4.1, 2.1))

i=1
for (i in 1:length(all_tables_list)){ #in normal circumstances
#for (i in 24:47){ #when plot has to be split, modify
  print(i)
  print(file_list[i])
  if (length(all_tables_list[[i]]) < 1) {
    barplot(0,0, main = substr(file_list[i],1, 13), ylim = c(0,1))
    text(0.15, 0.5, "no bars passed")
    text(0.1, 0.3, "the filtering threshold")
  }else{
    file = file_list[i]
    barplot(all_tables_list[[i]], #$file
            col = colors, 
            main = substr(file_list[i],1, 13), #file, 1, 13),
            xlab = "contig numbers", 
            ylab = "No. of genes", 
            cex.lab = 1.25)
  }
}

### LEGEND ###
#plot grid
par(mar = c(0, 0, 1, 0)) #bottom,left,top,right
plot(1,1, type = "n", bty="n",
     mar = c(1,1,1,1),
     main = "Legend", xpd = TRUE,
     xaxt = "n", yaxt = "n",
     xlab = c(" "), ylab = c(" "))

legend(0.65,1.4,  inset=c(-0.1, 0),
       legend = rev(order_vector), 
       fill = rev(colors))


#one plot
legend(ncol(all_tables_list[[i]])+0.25, 
       max(colSums(all_tables_list[[i]])), 
       legend = rev(order_vector), 
       fill = rev(colors))






