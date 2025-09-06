# GO_enrichment_barplot.R
# Purpose: Visualize top GO terms from g:Profiler results
# Author: Donovin Coles
# Date: [2019-09-23]

# Load libraries
library(ggplot2)
library(tidyverse)
library(viridisLite)
library(viridis)
library(gridExtra)

# Read g:Profiler results
gprofiler <- read.csv("Top20GOterms_gProfiler.csv")

# Prepare data
gprofiler <- gprofiler %>% 
  mutate(across(4:5, as.numeric)) %>%          # Ensure numeric columns
  mutate(source = as.factor(source)) %>%       # Ensure source is factor
  arrange(intersection_size, source) %>% 
  glimpse()

# Plot all GO terms
all_plot <- ggplot(gprofiler, aes(
    x = reorder(term_name, intersection_size),
    y = intersection_size,
    color = adjusted_p_value
  )) +
  geom_segment(aes(xend = term_name, yend = 0), color = "darkgrey", size = 1.3) +
  geom_point(size = 4) +
  facet_grid(source ~ ., scales = "free", drop = TRUE, space = "free") +
  scale_color_viridis(begin = 0.95, end = 0.2) +
  coord_flip() +
  theme_bw() +
  theme(legend.position = "top", axis.text = element_text(size = 13))

# Plot specific source (example: source "12")
gprofiler_12 <- gprofiler %>% filter(source == "12") %>% glimpse()

g12_plot <- ggplot(gprofiler_12, aes(
    x = reorder(term_name, intersection_size),
    y = intersection_size,
    color = adjusted_p_value
  )) +
  geom_segment(aes(xend = term_name, yend = 0), color = "darkgrey", size = 1.3) +
  geom_point(size = 4) +
  scale_color_viridis(begin = 0.95, end = 0.2) +
  ylim(c(0, 60)) +
  coord_flip() +
  theme_bw() +
  theme(legend.position = "none", axis.text = element_text(size = 13))

# Arrange plots together
grid.arrange(all_plot, g12_plot, nrow = 2)
