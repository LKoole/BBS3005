# Clear global environment
rm(list=ls())

# Set working directory
setwd("/Users/lisakoole/Desktop/Datasets")

#import data from file
library(tidyr)
file <- file.path(getwd(),"Pathway_statistics.txt")
dat <- read.delim(file, header = TRUE, sep = "\t", dec = ",", skip = 8)

# Separate column into two columns
dat2 <- dat %>% separate(X., c("X", "Percentage"), sep = "%")

# Convert "," into "." 
dat2[,5] <- as.numeric(gsub(",", ".", gsub("\\.", "", dat2$X)))

# Add column names
colnames(dat2) <- paste(c("Pathway", "positive.r", "measured.n", "total", "X", "Percentage", "Z.score", "p.value.perm"))

# Export data.frame to txt file
dat3 <- dat2[,c(1,2,3,4,5,7,8)]
write.table(dat3, "Pathway_statistics_formatted.txt", sep = "\t", dec = ".", quote = FALSE, row.names = FALSE, col.names = TRUE)
