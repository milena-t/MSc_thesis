#### use table from my own script to get cumulative pangenome size
library(ggplot2)
install.packages("fields")
library(fields)
install.packages("viridis")  
library("viridis") 

# home

setwd("/home/milena/megasync/masters_thesis/genome_sets_for_evaluation/long_sections/cluster1_results/bootstrap_curve")
setwd("/home/milena/megasync/masters_thesis/genome_sets_for_evaluation/single_STs/ST_50_results/bootstrap_curve")

# work      
setwd("/home/milena/megasync/masters_thesis/jejuni_1000_random_genomes/jejuni_1000_-i50_analysis/roary_outfiles")

cumulative_gene_content <- read.table("cumulative_gene_content_2.tsv", header = T, sep = '\t')

pangenome_no <- c(1:length(cumulative_gene_content$genome))
differences <- cumulative_gene_content$genomes_added

##################### MLST CLUSTERS COLORS ########################


## color by phylogenetic subgroup (20 clusters)
setwd("/media/milena/One Touch/mlst/jejuni1000_analysis")
clusters_strings <- scan("all_clusters_R.csv", what = "", sep = "\n")
# scan because the number of elements is different in each row 
clusters <- list()
len <- length(clusters)
# make list with strings of each cluster
for (i in 1:length(clusters_strings)){
  clusters[i] <- strsplit(clusters_strings[i], ',')
}
# names from the pangenome output (correct order)
setwd("/media/milena/One Touch/jejuni_1000_-i90")
genome_names <- (read.table("gene_presence_absence_genomes_reshuffled.Rtab"))[1,2:1000]
# find which cluster has each genome name and assign a color
cols <- rep("#ABABAB", length(genome_names)) #base color grey to see if any genomes are skipped
color_palette <- viridis(length(clusters_strings))

for (g in 1:length(genome_names)){
  #extract character string and remove file extension
  genome <- genome_names[[g]]
  genome <-sub('.1.gbff*', "", genome)
  for (i in 1:length(clusters_strings)){
    genome_in_cluster <- grepl(genome, clusters_strings[i], fixed = T)
    #print(genome_in_cluster)
    if (genome_in_cluster){
      #print(color_palette[i])
      cols[g] <- color_palette[i]
      break # go to next genome, do not continue with checking the following cluster_strings 
    }
  }
}
#find the first instance of each color in the color vector
x_first_color = numeric(length(color_palette))
for (i in 1:length(color_palette)){
  x_first_color[i]<- match(color_palette[i], cols)
}


##################### PLOT DIFFERENCES CURVE ######################

# if colors not by phylogenetic subgroup
diff_q <- quantile(differences[75:length(differences)], 0.975, na.rm = T)
cols <- c(rep("cornflowerblue", length(pangenome_no)))
cols[differences>diff_q] <- "firebrick"
cols[1:75] <- "darkgrey"

plot(1:length(differences), differences, 
     ylim = c(0, 150),
     xlab = "no. genomes", 
     ylab = "no. of new gene families added", 
     main = "No. of Genes added for each genome ST 50",
     pch = 19, col = cols)
#yline(diff_q, lty = "dashed")
#yline(median(differences[500:999]), lty = "dashed")
#text(35, 350, "first 100 genomes", col = "darkgrey", cex = 0.8)
text(500, diff_q+100, "97.5% quantile", col = "firebrick")#, cex = -50)
# 
#### fit curve ####  
fit_function <- function(x, a, b){
  return((1*a)/x+b)
}
curve(fit_function(x, 1000, median(cumulative_gene_content$genomes_added)), add = T)

### if colors by cluster:
#plot first instance of new color when colored by phylogenetic cluster
for (i in 1:length(x_first_color)){
  xline(x_first_color[i], col = color_palette[i], lty = 'dashed')
}


##################### PLOT PANGENOME CURVE #########################


plot(pangenome_no, cumulative_gene_content$pangenome_size, type = "l", lwd = 2, 
     ylim = c(0, 15000), 
     xlab = "No. of genomes", ylab = "No. of gene families",
     main = "Pan-genome curves", 
)
lines(pangenome_no, cumulative_gene_content1$pangenome_size, lwd = 2)
# mark jumps in gene content I found
for (i in indices){
  xline(i, col = "darkgrey", lty = "dashed")
}

##################### PLOT BOOTSTRAPPED PANGENOME #####################

plot(0,0, type = "n", 
     ylim = c(0, 13000), 
     xlim = c(0, 1000),
     xlab = "No. of genomes in the pan-genome", ylab = "No. of gene families",
     main = "Pan-genome curves (with random permutations)", 
     col = "seagreen", cex.main = 1.5, cex.lab = 1.4
)

##### ST48 #####
if(T){
setwd("/home/milena/megasync/masters_thesis/genome_sets_for_evaluation/single_STs/ST_48_results/bootstrap_curve")
cumulative_gene_content1 <- read.table("cumulative_gene_content_1.tsv", header = T, sep = '\t')
cumulative_gene_content2 <- read.table("cumulative_gene_content_2.tsv", header = T, sep = '\t')
cumulative_gene_content3 <- read.table("cumulative_gene_content_3.tsv", header = T, sep = '\t')
cumulative_gene_content4 <- read.table("cumulative_gene_content_4.tsv", header = T, sep = '\t')
cumulative_gene_content5 <- read.table("cumulative_gene_content_5.tsv", header = T, sep = '\t')
cumulative_gene_content6 <- read.table("cumulative_gene_content_6.tsv", header = T, sep = '\t')
cumulative_gene_content7 <- read.table("cumulative_gene_content_7.tsv", header = T, sep = '\t')
cumulative_gene_content8 <- read.table("cumulative_gene_content_8.tsv", header = T, sep = '\t')
cumulative_gene_content9 <- read.table("cumulative_gene_content_9.tsv", header = T, sep = '\t')
cumulative_gene_content <- read.table("cumulative_gene_content.tsv", header = T, sep = '\t')
pangenome_no <- c(1:length(cumulative_gene_content$genome))

# make mean of all curves
cumulative_means2 <- cumulative_gene_content$pangenome_size
for (i in 1:length(cumulative_gene_content$pangenome_size)){
  nums <- c(
    cumulative_gene_content$pangenome_size[i],
    cumulative_gene_content1$pangenome_size[i],
    cumulative_gene_content2$pangenome_size[i],
    cumulative_gene_content3$pangenome_size[i],
    cumulative_gene_content4$pangenome_size[i],
    cumulative_gene_content5$pangenome_size[i],
    cumulative_gene_content6$pangenome_size[i],
    cumulative_gene_content7$pangenome_size[i],
    cumulative_gene_content8$pangenome_size[i],
    cumulative_gene_content9$pangenome_size[i]
  )
  cumulative_means2[i]<- mean(nums)
}

lines(pangenome_no, cumulative_gene_content1$pangenome_size, lwd = 2, col = "darkseagreen3")
lines(pangenome_no, cumulative_gene_content2$pangenome_size, lwd = 2, col = "darkseagreen3")
lines(pangenome_no, cumulative_gene_content3$pangenome_size, lwd = 2, col = "darkseagreen3")
lines(pangenome_no, cumulative_gene_content4$pangenome_size, lwd = 2, col = "darkseagreen3")
lines(pangenome_no, cumulative_gene_content5$pangenome_size, lwd = 2, col = "darkseagreen3")
lines(pangenome_no, cumulative_gene_content6$pangenome_size, lwd = 2, col = "darkseagreen3")
lines(pangenome_no, cumulative_gene_content7$pangenome_size, lwd = 2, col = "darkseagreen3")
lines(pangenome_no, cumulative_gene_content8$pangenome_size, lwd = 2, col = "darkseagreen3")
lines(pangenome_no, cumulative_gene_content9$pangenome_size, lwd = 2, col = "darkseagreen3")
lines(pangenome_no, cumulative_gene_content$pangenome_size, lwd = 2, col = "darkseagreen3")
lines(pangenome_no, cumulative_means2, lwd = 2, col = "seagreen")
}

##### ST50 #####
if(T){
setwd("/home/milena/megasync/masters_thesis/genome_sets_for_evaluation/single_STs/ST_50_results/bootstrap_curve")
cumulative_gene_content1 <- read.table("cumulative_gene_content_1.tsv", header = T, sep = '\t')
cumulative_gene_content2 <- read.table("cumulative_gene_content_2.tsv", header = T, sep = '\t')
cumulative_gene_content3 <- read.table("cumulative_gene_content_3.tsv", header = T, sep = '\t')
cumulative_gene_content4 <- read.table("cumulative_gene_content_4.tsv", header = T, sep = '\t')
cumulative_gene_content5 <- read.table("cumulative_gene_content_5.tsv", header = T, sep = '\t')
cumulative_gene_content6 <- read.table("cumulative_gene_content_6.tsv", header = T, sep = '\t')
cumulative_gene_content7 <- read.table("cumulative_gene_content_7.tsv", header = T, sep = '\t')
cumulative_gene_content8 <- read.table("cumulative_gene_content_8.tsv", header = T, sep = '\t')
cumulative_gene_content9 <- read.table("cumulative_gene_content_9.tsv", header = T, sep = '\t')
cumulative_gene_content <- read.table("cumulative_gene_content.tsv", header = T, sep = '\t')
pangenome_no <- c(1:length(cumulative_gene_content$genome))

# make mean of all curves
cumulative_means1 <- cumulative_gene_content$pangenome_size
for (i in 1:length(cumulative_gene_content$pangenome_size)){
  nums <- c(
    cumulative_gene_content$pangenome_size[i],
    cumulative_gene_content1$pangenome_size[i],
    cumulative_gene_content2$pangenome_size[i],
    cumulative_gene_content3$pangenome_size[i],
    cumulative_gene_content4$pangenome_size[i],
    cumulative_gene_content5$pangenome_size[i],
    cumulative_gene_content6$pangenome_size[i],
    cumulative_gene_content7$pangenome_size[i],
    cumulative_gene_content8$pangenome_size[i],
    cumulative_gene_content9$pangenome_size[i]
  )
  cumulative_means1[i]<- mean(nums)
}

lines(pangenome_no, cumulative_gene_content1$pangenome_size, lwd = 2, col = "thistle3")
lines(pangenome_no, cumulative_gene_content2$pangenome_size, lwd = 2, col = "thistle3")
lines(pangenome_no, cumulative_gene_content3$pangenome_size, lwd = 2, col = "thistle3")
lines(pangenome_no, cumulative_gene_content4$pangenome_size, lwd = 2, col = "thistle3")
lines(pangenome_no, cumulative_gene_content5$pangenome_size, lwd = 2, col = "thistle3")
lines(pangenome_no, cumulative_gene_content6$pangenome_size, lwd = 2, col = "thistle3")
lines(pangenome_no, cumulative_gene_content7$pangenome_size, lwd = 2, col = "thistle3")
lines(pangenome_no, cumulative_gene_content8$pangenome_size, lwd = 2, col = "thistle3")
lines(pangenome_no, cumulative_gene_content9$pangenome_size, lwd = 2, col = "thistle3")
lines(pangenome_no, cumulative_gene_content$pangenome_size, lwd = 2, col = "thistle3")
lines(pangenome_no, cumulative_means1, lwd = 2, col = "slateblue")
}

##### Group 1 #####
if(T){
setwd("/home/milena/megasync/masters_thesis/genome_sets_for_evaluation/long_sections/cluster1_results/bootstrap_curve")
cumulative_gene_content1 <- read.table("cumulative_gene_content_1.tsv", header = T, sep = '\t')
cumulative_gene_content2 <- read.table("cumulative_gene_content_2.tsv", header = T, sep = '\t')
cumulative_gene_content3 <- read.table("cumulative_gene_content_3.tsv", header = T, sep = '\t')
cumulative_gene_content4 <- read.table("cumulative_gene_content_4.tsv", header = T, sep = '\t')
cumulative_gene_content5 <- read.table("cumulative_gene_content_5.tsv", header = T, sep = '\t')
cumulative_gene_content6 <- read.table("cumulative_gene_content_6.tsv", header = T, sep = '\t')
cumulative_gene_content7 <- read.table("cumulative_gene_content_7.tsv", header = T, sep = '\t')
cumulative_gene_content8 <- read.table("cumulative_gene_content_8.tsv", header = T, sep = '\t')
cumulative_gene_content9 <- read.table("cumulative_gene_content_9.tsv", header = T, sep = '\t')
cumulative_gene_content <- read.table("cumulative_gene_content.tsv", header = T, sep = '\t')
pangenome_no <- c(1:length(cumulative_gene_content$genome))

# make mean of all curves
cumulative_means2 <- cumulative_gene_content$pangenome_size
for (i in 1:length(cumulative_gene_content$pangenome_size)){
  nums <- c(
    cumulative_gene_content$pangenome_size[i],
    cumulative_gene_content1$pangenome_size[i],
    cumulative_gene_content2$pangenome_size[i],
    cumulative_gene_content3$pangenome_size[i],
    cumulative_gene_content4$pangenome_size[i],
    cumulative_gene_content5$pangenome_size[i],
    cumulative_gene_content6$pangenome_size[i],
    cumulative_gene_content7$pangenome_size[i],
    cumulative_gene_content8$pangenome_size[i],
    cumulative_gene_content9$pangenome_size[i]
  )
  cumulative_means2[i]<- mean(nums)
}


lines(pangenome_no, cumulative_gene_content1$pangenome_size, lwd = 2, col = "darkslategray3")
lines(pangenome_no, cumulative_gene_content2$pangenome_size, lwd = 2, col = "darkslategray3")
lines(pangenome_no, cumulative_gene_content3$pangenome_size, lwd = 2, col = "darkslategray3")
lines(pangenome_no, cumulative_gene_content4$pangenome_size, lwd = 2, col = "darkslategray3")
lines(pangenome_no, cumulative_gene_content5$pangenome_size, lwd = 2, col = "darkslategray3")
lines(pangenome_no, cumulative_gene_content6$pangenome_size, lwd = 2, col = "darkslategray3")
lines(pangenome_no, cumulative_gene_content7$pangenome_size, lwd = 2, col = "darkslategray3")
lines(pangenome_no, cumulative_gene_content8$pangenome_size, lwd = 2, col = "darkslategray3")
lines(pangenome_no, cumulative_gene_content9$pangenome_size, lwd = 2, col = "darkslategray3")
lines(pangenome_no, cumulative_gene_content$pangenome_size, lwd = 2, col = "darkslategray3")
lines(pangenome_no, cumulative_means2, lwd = 2, col = "darkslategray4")
}

##### Group 2 #####
if(T){
  setwd("/home/milena/megasync/masters_thesis/genome_sets_for_evaluation/long_sections/cluster2_results/bootstrap_curve")
  cumulative_gene_content1 <- read.table("cumulative_gene_content_1.tsv", header = T, sep = '\t')
  cumulative_gene_content2 <- read.table("cumulative_gene_content_2.tsv", header = T, sep = '\t')
  cumulative_gene_content3 <- read.table("cumulative_gene_content_3.tsv", header = T, sep = '\t')
  cumulative_gene_content4 <- read.table("cumulative_gene_content_4.tsv", header = T, sep = '\t')
  cumulative_gene_content5 <- read.table("cumulative_gene_content_5.tsv", header = T, sep = '\t')
  cumulative_gene_content6 <- read.table("cumulative_gene_content_6.tsv", header = T, sep = '\t')
  cumulative_gene_content7 <- read.table("cumulative_gene_content_7.tsv", header = T, sep = '\t')
  cumulative_gene_content8 <- read.table("cumulative_gene_content_8.tsv", header = T, sep = '\t')
  cumulative_gene_content9 <- read.table("cumulative_gene_content_9.tsv", header = T, sep = '\t')
  cumulative_gene_content <- read.table("cumulative_gene_content.tsv", header = T, sep = '\t')
  pangenome_no <- c(1:length(cumulative_gene_content$genome))
  
  # make mean of all curves
  cumulative_means2 <- cumulative_gene_content$pangenome_size
  for (i in 1:length(cumulative_gene_content$pangenome_size)){
    nums <- c(
      cumulative_gene_content$pangenome_size[i],
      cumulative_gene_content1$pangenome_size[i],
      cumulative_gene_content2$pangenome_size[i],
      cumulative_gene_content3$pangenome_size[i],
      cumulative_gene_content4$pangenome_size[i],
      cumulative_gene_content5$pangenome_size[i],
      cumulative_gene_content6$pangenome_size[i],
      cumulative_gene_content7$pangenome_size[i],
      cumulative_gene_content8$pangenome_size[i],
      cumulative_gene_content9$pangenome_size[i]
    )
    cumulative_means2[i]<- mean(nums)
  }

  
  lines(pangenome_no, cumulative_gene_content1$pangenome_size, lwd = 2, col = "goldenrod")
  lines(pangenome_no, cumulative_gene_content2$pangenome_size, lwd = 2, col = "goldenrod")
  lines(pangenome_no, cumulative_gene_content3$pangenome_size, lwd = 2, col = "goldenrod")
  lines(pangenome_no, cumulative_gene_content4$pangenome_size, lwd = 2, col = "goldenrod")
  lines(pangenome_no, cumulative_gene_content5$pangenome_size, lwd = 2, col = "goldenrod")
  lines(pangenome_no, cumulative_gene_content6$pangenome_size, lwd = 2, col = "goldenrod")
  lines(pangenome_no, cumulative_gene_content7$pangenome_size, lwd = 2, col = "goldenrod")
  lines(pangenome_no, cumulative_gene_content8$pangenome_size, lwd = 2, col = "goldenrod")
  lines(pangenome_no, cumulative_gene_content9$pangenome_size, lwd = 2, col = "goldenrod")
  lines(pangenome_no, cumulative_gene_content$pangenome_size, lwd = 2, col = "goldenrod")
  lines(pangenome_no, cumulative_means2, lwd = 2, col = "goldenrod4")
}

# plot comparative randomly selected line
setwd("/home/milena/megasync/masters_thesis/jejuni_1000_random_genomes/jejuni_1000_-i50_analysis/roary_outfiles")
cumulative_gene_content <- read.table("cumulative_gene_content.tsv", header = T, sep = '\t')
pangenome_no <- c(1:length(cumulative_gene_content$genome))
lines(pangenome_no, cumulative_gene_content$pangenome_size, lwd = 2, col = "gray")

legend(x = 'bottomright', legend = c("group 1", "group 2", "ST 50", "ST 48", "randomly selected genomes"), 
      col = c( "darkslategray3", "goldenrod", "thistle3", "darkseagreen3" ,"gray"), lwd = 2, cex = 1.25)


################### FIT POWER LAW ##################
# the power law fits to the no of genes added, so the differences curve.
# I will fit it to the mean of the permutations in the differences curve

## Group 1

cumulative_means1 <- cumulative_gene_content$pangenome_size
for (i in 1:length(cumulative_gene_content$pangenome_size)){
  nums <- c(
    cumulative_gene_content$genomes_added[i],
    cumulative_gene_content1$genomes_added[i],
    cumulative_gene_content2$genomes_added[i],
    cumulative_gene_content3$genomes_added[i],
    cumulative_gene_content4$genomes_added[i],
    cumulative_gene_content5$genomes_added[i],
    cumulative_gene_content6$genomes_added[i],
    cumulative_gene_content7$genomes_added[i],
    cumulative_gene_content8$genomes_added[i],
    cumulative_gene_content9$genomes_added[i]
  )
  cumulative_means1[i]<- mean(nums)
}

no_genomes <- seq(1, length(cumulative_means1))
model1 <- lm(cumulative_means1~log(no_genomes))
model1 <- lm(log(cumulative_means1)~no_genomes)
model1 <- lm(log(cumulative_means1)~log(no_genomes))
summary(model1)
model1_exponential <- exp(predict(model1, list(no_genomes = no_genomes)))
plot(no_genomes, cumulative_means1 ,pch = 16, log = "xy",
     #ylim = c(0, 2000), 
     xlab = "avg. no. of new gene families added", ylab = "No. of gene families",
     main = "Group 1 power-law fit", 
     col = "slateblue"
)
lines(no_genomes, model1_exponential, lwd = 2)

fit_function <- function(x, a, k){
  return(k+x^(-a))
}
for (i in 0:1000){
  curve(fit_function(x, i, median(cumulative_means1)), add = T)
}
curve(fit_function(x, 1, median(cumulative_means1)), add = T)

### if colors by cluster:
#plot first instance of new color when colored by phylogenetic cluster
for (i in 1:length(x_first_color)){
  xline(x_first_color[i], col = color_palette[i], lty = 'dashed')
}

################### PLOT MULTIPLE PANGENOME CURVES IN ONE PLOT ###############

setwd("/home/milena/mega/masters_thesis/genome_sets_for_evaluation/long_sections/cluster1_results")
cumulative_gene_content1 <- read.table("cumulative_gene_content.tsv", header = T, sep = '\t')
pangenome_no1 <- c(1:length(cumulative_gene_content1$genome))
differences1 <- cumulative_gene_content1$genomes_added
setwd("/home/milena/mega/masters_thesis/genome_sets_for_evaluation/long_sections/cluster2_results")
cumulative_gene_content2 <- read.table("cumulative_gene_content.tsv", header = T, sep = '\t')
pangenome_no2 <- c(1:length(cumulative_gene_content2$genome))
differences2 <- cumulative_gene_content2$genomes_added

plot(0, 0, ylim = c(0, 12000), xlim = c(0, 600), type = 'n',
     xlab = "No. of genomes", ylab = "No. of gene families",
     main = "Pan-genome curves for two C. jejuni datasets")

lines(pangenome_no1, cumulative_gene_content1$pangenome_size, type = "l", lwd = 2, 
        col = "seagreen")
lines(pangenome_no2, cumulative_gene_content2$pangenome_size, type = "l", lwd = 2, 
      col = "slateblue")

legend(x = 'bottomright', legend = c("Group 1", "Group 2"), col = c("seagreen", "slateblue"), lwd = 2)



################ GET TOO-HIGH GENOME NAMES ####################

too_high_addition_index <- cols=="firebrick"
too_high_genome_names <- cumulative_gene_content$genome[too_high_addition_index]
# remove file ending
for (i in 1:length(too_high_genome_names)){
  too_high_genome_names[i] <- sub('.1.gbff*', "", too_high_genome_names[i])
}
print(too_high_genome_names)

#### get nomral genomes as negative control

x_cos_blue <- pangenome_no[cols=="cornflowerblue"]
normal_addition_index <- sample(x_cos_blue, length(too_high_genome_names))
normal_genome_names <- cumulative_gene_content$genome[normal_addition_index]
print(normal_genome_names)
















