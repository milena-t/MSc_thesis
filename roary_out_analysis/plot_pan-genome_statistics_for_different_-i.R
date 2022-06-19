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
barplot(matrix)

axis(1, at = seq_identity_int)
