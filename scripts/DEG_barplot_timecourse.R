# DEG_barplot_timecourse.R
# Purpose: Visualize the number of differentially expressed genes (DEGs) over a time course
# Author: Donovin Coles

# Load libraries
library(ggplot2)

# Read data
DEG <- read.csv("SDEG.csv")
head(DEG)

# Create bar plot of DEGs over time
p <- ggplot(DEG, aes(x = HPI, y = DEG, fill = Regulation)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(
    x = "\nTime after inoculation (h)", 
    y = "Differentially expressed genes \n"
  ) +
  scale_y_continuous(limits = c(0, 10000)) +
  theme_minimal() +
  scale_fill_grey() +
  theme(
    text = element_text(size = 16)
  )

# Print plot
print(p)
