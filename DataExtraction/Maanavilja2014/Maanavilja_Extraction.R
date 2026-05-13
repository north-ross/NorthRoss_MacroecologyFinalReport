mosses <- read.csv("PCtCoverMosses.csv", header = FALSE)
shrubHedge <- read.csv("PctCoverShrubHedge.csv", header = FALSE)

# Convert SD to a field
library(dplyr)
library(tidyr)
# Pivot
mosses <- mosses %>%
  pivot_wider(
    names_from = 4,
    values_from = 2
  )

mosses$n <- rep(9, nrow(mosses))
mosses$mean = mosses$` Mean`
# Convert SD to actual amount by subtracting mean
mosses[,5] = mosses[,5] - mosses[,4]

# SD to var
mosses$var = mosses[,5] ^2

# Sum of means and variances
mosses %>%
  group_by(V5) %>%
  summarise(meansum = sum(mean), 
            varsum_SE = sqrt(sum(var)/mean(n)), 
            n = mean(n))

#####################
# Repeat for others

# Pivot
shrubHedge <- shrubHedge %>%
  pivot_wider(
    names_from = 4,
    values_from = 2
  )

shrubHedge$n <- rep(9, nrow(shrubHedge))
shrubHedge$mean = shrubHedge$` Mean`
# Convert SD to actual amount by subtracting mean
shrubHedge$` SD` = shrubHedge$` SD` - shrubHedge$mean

# SD to var
shrubHedge$var = shrubHedge$` SD` ^2

# Sum of means and variances
shrubHedge %>%
  group_by(V5, V6) %>%
  summarise(meansum = sum(mean), 
            varsum_SE = sqrt(sum(var)/mean(n)), 
            n = mean(n))

