setwd(" folder location")

#input data
name1 <- read.csv ("file name.csv")

name2 <- read.csv ("file name.csv")

install.packages("dplyr")
library("dplyr")


#merge
merged <- left_join(name1, name2, by="gene_id")


#write files
write.csv(merged, "output file name.csv")
