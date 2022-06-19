# make pangenome curve for different blastp thresholds in roary
library(fields)
library("viridis")
library(RColorBrewer)

setwd('/home/milena/megasync/masters_thesis/jejuni_1000_random_genomes/jejuni_1000_different_blastp_thresholds')
cumulative_pangneome_size_files = c(
  'roary_results_-i40_1645039083/cumulative_gene_content.tsv',
  'roary_results_-i45_1645036244/cumulative_gene_content.tsv',
  'roary_results_-i50_1645033345/cumulative_gene_content.tsv',
  'roary_results_-i55_1645030419/cumulative_gene_content.tsv',
  'roary_results_-i60_1645027528/cumulative_gene_content.tsv',
  'roary_results_-i65_1645024647/cumulative_gene_content.tsv',
  'roary_results_-i70_1645021689/cumulative_gene_content.tsv',
  'roary_results_-i75_1645018784/cumulative_gene_content.tsv',
  'roary_results_-i80_1645015900/cumulative_gene_content.tsv',
  'roary_results_-i85_1645013020/cumulative_gene_content.tsv',
  'roary_results_-i90_1645010179/cumulative_gene_content.tsv',
  'roary_results_-i95_1645007285/cumulative_gene_content.tsv')

seq_identity <- c('40%', '45%', '50%', '55%', '60%', '65%', '70%', '75%', '80%', '85%', '90%', '95%')

color_palette <- viridis((length(cumulative_pangneome_size_files)))

############## PLOT PANGENOME CURVE ################

par(mar=c(5, 4, 4, 8), xpd=TRUE) #make space for the legend to the right of the plot
par(mar=c(5, 4, 4, 3), xpd=TRUE)
# make empty plot with correct dimensions
plot(0, 0, ylim = c(0, 20000), xlim = c(0, 1000), type = 'n',
     xlab = "No. of genomes", ylab = "No. of gene families",
     main = "Pan-genome curves for different BLASTP thresholds")

indices_to_plot <- c(3,11) #only 90% and 50%
indices_to_plot <- rev(1:length(cumulative_pangneome_size_files)) #all
# add lines iteratively
for (i in indices_to_plot){
  #print(i)
  cumulative_gene_content <- read.table(cumulative_pangneome_size_files[i], header = T, sep = '\t')
  pangenome_no <- c(1:length(cumulative_gene_content$genome))
  
 lines(pangenome_no, cumulative_gene_content$pangenome_size, type = "l", lwd = 2, 
       col = color_palette[i])
}
legend(x = 'topright', inset=c(-0.3, 0), legend = rev(seq_identity), col = rev(color_palette), lwd = 2)
legend(x = 'bottomright', legend = rev(seq_identity[indices_to_plot]), col = rev(color_palette[indices_to_plot]), lwd = 2)


################ PLOT PANGENOME SIZE ################

pangenome_sizes <- c()
seq_identity_int <- seq(40, 95, 5)
for (i in 1:length(cumulative_pangneome_size_files)){
  #print(i)
  cumulative_gene_content <- read.table(cumulative_pangneome_size_files[i], header = T, sep = '\t')
  pangenome_no <- c(1:length(cumulative_gene_content$genome))
  print(cumulative_gene_content$pangenome_size[length(cumulative_gene_content$pangenome_size)])
  pangenome_sizes <- append(pangenome_sizes, cumulative_gene_content$pangenome_size[length(cumulative_gene_content$pangenome_size)])
}

plot(0, 0, ylim = c(10000, 20000), xlim = c(44, 96), xaxt = 'n', type = 'n',
     xlab = "% Seq. Identity", ylab = "No. of gene families",
     main = "Pan-genome size for different BLASTP thresholds")
axis(1, at = seq_identity_int)

points(seq_identity_int, pangenome_sizes, pch = 19)#, col = color_palette)


################### pan-genome composition ###################

# create bar plot to show different pan-genome compositions for different -i thresholds
# based on the roary outfiles "summary_statistics.txt

# order in the vectors:
# core genes (appear in at least 99\% of genomes), 
# and accessory genes split in three different categories: 
#   Soft core genes (appear in 95\% to 99\% of genomes), 
#   shell genes (appear in 15\% to 95\%), 
#   and cloud genes (appear in 0\% to 15\%)
# pan-genome size at the end

normalize <- function(i90){
  i90_normalized <- i90/i90[5]
  return(i90_normalized[1:4])
}

i95 <- c(738, 327, 1087, 15661, 17813)
i90 <- c(774, 307, 1036, 14161, 16278)
i85 <- c(798, 291, 1036, 13072, 15197)
i80 <- c(840, 252, 1024, 12394, 14510)
i75 <- c(827, 259, 1049, 11801, 13936)
i70 <- c(843, 237, 1060, 11264, 13404)
i65 <- c(851, 233, 1063, 11139, 13286)
i60 <- c(860, 220, 1048, 10821, 12949)
i55 <- c(843, 231, 1087, 10749, 12910)
i50 <- c(831, 226, 1106, 10780, 12943)
i45 <- c(848, 242, 1025, 10596, 12711)
i40 <- c(844, 219, 1101, 10629, 12793)

matrix <- cbind(normalize(i45), normalize(i50), normalize(i55), 
                normalize(i60), normalize(i65), normalize(i70), 
                normalize(i75), normalize(i80), normalize(i85),
                normalize(i90), normalize(i95))

cols = rev(brewer.pal(4, 'Blues'))

barplot(matrix, main = "Pan-genome composition for different clustering thresholds", 
        names.arg = c("45%", "50%", "55%", "60%", "65%", "70%", "75%", "80%", "85%", "90%", "95%"),
        ylab = "proportion of the pan-genome", xlab = "sequence identity threshold", 
        col = cols)

legend("topright", legend = rev(c("core", "soft core", "shell", "cloud")), fill = rev(cols))

# pairwise comparison of 45% and 90%

t.test_45v90 <- t.test(normalize(i45), normalize(i90), alternative = "two.sided", paired = T)
t.test_45v90
 

################ PLOT MAIN DATASETS ################

group1 <- c(890, 310, 900, 7260, 9360)
group2 <- c(687, 380, 932, 8831, 10830)
st48 <- c(767, 489, 580, 4007, 5843)
st50 <- c(777, 440, 758, 4314, 6289)

matrix <- cbind(normalize(group1), normalize(group2),
                normalize(st50), normalize(st48))
matrix <- cbind(normalize(st50), normalize(st48))

cols = rev(brewer.pal(4, 'Blues'))
cols <- rev(c("#b3cac7", "#82ada8", "#578e87", "#3a6c65"))

barplot(matrix, main = "Pan-genome composition", 
        names.arg = c("Group 1", "Group 2", "ST 50", "ST 48"),
        ylab = "proportion of the pan-genome", xlab = "dataset", 
        col = cols, cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.25, cex.names = 1.25)

barplot(matrix, main = "Pan-genome composition", 
        names.arg = c("ST 50", "ST 48"),
        ylab = "proportion of the pan-genome", xlab = "dataset", 
        col = cols, cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.25, cex.names = 1.25)

legend("topright", legend = rev(c("core", "soft core", "shell", "cloud")), fill = rev(cols), cex = 1.5)







