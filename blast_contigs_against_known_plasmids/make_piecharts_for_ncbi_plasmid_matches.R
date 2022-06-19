# make pie charts for ncbi plasmid comparisons

datasets <- c("Group 1", "Group 2", "ST 48", "ST 50")
no_contigs <- c(77, 135, 147, 431)
no_not_matched <- c(2, 4, 23, 213)
no_matched <- no_contigs - no_not_matched
cols <- c("#b3cac7", "#50938a")

par(mfrow = c(1,4))

#################### PIE CHART ################

for (i in 1:4){
  pie(c(no_matched[i], no_not_matched[i]), 
      labels = c(as.character(no_matched[i]), as.character(no_not_matched[i])), 
      col = cols, main = datasets[i], cex.main = 2.5,
      sub = paste(as.character(no_contigs[i]), "contigs"), cex.sub = 2,
      radius = 1, cex=1.5)
  if (i == 2){
    legend("bottomright", c("have matches", "have no matches"), cex = 2, fill = cols)
    
  }
}

