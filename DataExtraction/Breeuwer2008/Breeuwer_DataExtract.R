
# water table depths - from web plot digitizer
CWT <- read.csv('CWT.csv', header=FALSE)
FWT <- read.csv('FWT.csv', header=FALSE)

mean(CWT[,2])
length(CWT[,2])
sd(CWT[,2])

mean(FWT[,2])
length(FWT[,2])
sd(FWT[,2])


# Moss growth rates from Table 4
moss <- data.frame(
  Species = c(1,2,3,1,2,3),
  # control vs treatment
  Group = c("C", "C", "C", "T", "T", "T"),
  means = c(640.7, 312.1, 47.0, 113.2, 645.1, 47.4),
  SE = c(121.4, 49.5, 15.3, 52.1, 99.2, 11.5),
  n = c(6,6,6,9,9,9)
)

# Moss cover change rates from Table 4
moss <- data.frame(
  Species = c(1,2,3,1,2,3),
  # control vs treatment
  Group = c("C", "C", "C", "T", "T", "T"),
  means = c(-20.2, 16.3, 3.5, -10.2, 11.1, -3.1),
  SE = c(7,7.4,2.7,3.7,3,2.5),
  n = c(6,6,6,9,9,9)
)
# Convert SE to Var
moss$var = moss$SE^2 * moss$n

# Report summed variance and means by control group
library(dplyr)
options(pillar.sigfig = 5)
moss %>%
  group_by(Group) %>%
  summarise(meansum = sum(means), varsum = sqrt(sum(var)/mean(n)), n = mean(n))


